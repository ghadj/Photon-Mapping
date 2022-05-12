#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "libs/SDLauxiliary.h"
#include "TestModel.h"

#include "photon_tracer.h"
#include "shapes.h"
#include "kdtree.h"
#include "utils.h"

using namespace std;

using glm::vec3;
using glm::mat3;
using glm::ivec2;
using glm::vec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 100;
const int SCREEN_HEIGHT = 100;
SDL_Surface* screen;
vector<Triangle> triangles; // all the trianlges of the scene
vector<Sphere> spheres;
int t; // time

// Camera
const float focalLength = SCREEN_HEIGHT;
vec3 cameraPos(0, 0, -3); // position of the camera

// Light
const vec3 lightPos(0, -0.5, -0.4);
const vec3 lightColor = 5.0f * vec3(1, 1, 0.95);
const vec3 lightPower = 400.f * vec3(1, 1, 1);

const int k = 500; // nearest photons
const float filter_const = 1.2;


const int nPhotons = 50000;
vector<Photon> photons;
KdTree photonmap;

// ----------------------------------------------------------------------------
// FUNCTIONS

vec3 DirectLight(const Intersection&);
void RayTrace();
void emitPhotons(int nPhotons);
void trace_photon(Photon& p);
void DrawPhoton(Photon p);
void DrawPhotonPath(Photon p);
ivec2 VertexShader(vec3 p);

int main(int argc, char* argv[])
{
    int t2; float dt;
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    LoadTestModel(triangles, spheres);
    
    emitPhotons(nPhotons);
    photonmap.setPhotons(photons.data(), photons.size());
    photonmap.buildTree();
    
    t2 = SDL_GetTicks();
    dt = float(t2-t);
    t = t2;
    cout << "Photonmap creation time: " << dt << " ms." << endl;
    cout<<"Photon-map size: "<<photons.size()<<'\n';

    RayTrace();
    
    // Compute frame time:
    t2 = SDL_GetTicks();
    dt = float(t2-t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;
    Uint8* keystate = SDL_GetKeyState( 0 );

    SDL_SaveBMP( screen, "screenshot.bmp" );
	
	return 0;
}

vec3 getRadianceEstimate(Intersection i)
{
    // Jensen eq. (11)
    vec3 color = vec3(0,0,0);
    vec3 delta_phi;
    float wpc; // weight
    float dp;  // distance between photon and point
    float r_sqr = 0.0; // max sqr distance from k-th photon

    // Get k-nearest photons
    vector<NeighborPhoton> neighbor_photons = photonmap.searchKNearest(i.position, 
                                                                       k, r_sqr);
    for(int p=0; p<neighbor_photons.size(); p++)
    {
        // Cone-filter [Jensen96c]
        dp = neighbor_photons[p].dist;
        wpc = 1 - dp/(k*sqrt(r_sqr));
        
        // Photon power
        delta_phi = max(glm::dot(-photons[neighbor_photons[p].index].getNormal(), 
                        triangles[i.triangleIndex].normal), 0.0f) *
                        photons[neighbor_photons[p].index].getEnergy();  
        color += wpc*delta_phi;
    }
    color /= (1 - 2/(3*filter_const) * PI * r_sqr);

    color += DirectLight(i); // Direct light
    color *= triangles[i.triangleIndex].color;

    return color;
}

void RayTrace()
{
	if(SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

    SDL_FillRect(screen, 0, 0);
    Intersection closestIntersection, intersection_refract, intersection_reflect;

	for(int y=0; y<SCREEN_HEIGHT; ++y)
		for(int x=0; x<SCREEN_WIDTH; ++x)
		{
            vec3 dir(x-SCREEN_WIDTH/2, y-SCREEN_HEIGHT/2, focalLength);
            if(ClosestIntersection(cameraPos, dir, triangles, spheres, closestIntersection))
            {
                if(closestIntersection.sphereIndex >= 0)
                {
                    refract(spheres[closestIntersection.sphereIndex],
                                  glm::normalize(closestIntersection.position-cameraPos), triangles, spheres,
                                  closestIntersection, intersection_refract);
                    
                    float f = fresnel(glm::normalize(closestIntersection.position-cameraPos),
                                      glm::normalize(closestIntersection.position - 
                                                     spheres[closestIntersection.sphereIndex].center));
                    float c = max(0.0, min(1.0, 1.5*f));

                    vec3 dir_re = reflect(glm::normalize(closestIntersection.position-cameraPos),
                                      glm::normalize(closestIntersection.position - 
                                                     spheres[closestIntersection.sphereIndex].center));
                    ClosestIntersection(closestIntersection.position, dir_re, triangles, spheres, intersection_reflect);

                    PutPixelSDL(screen, x, y, (1-c)*getRadianceEstimate(intersection_refract)+
                                            (c)*getRadianceEstimate(intersection_reflect));
                }
                else
                    PutPixelSDL(screen, x, y, getRadianceEstimate(closestIntersection));
            }
		}

	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

vec3 DirectLight(const Intersection& i)
{
    // Distance of object from light source
    float r = glm::distance(lightPos, i.position); 

    // Direction from surface to the light source
    vec3 r_hat = glm::normalize(lightPos - i.position); 

    // Normal pointing out from the surface
    vec3 n_hat = triangles[i.triangleIndex].normal; 
    
    Intersection j;
    if(ClosestIntersection(lightPos, glm::normalize(i.position - lightPos), 
                           triangles, spheres, j))
    {
        if(j.sphereIndex >= 0 ||
           (glm::distance(lightPos, j.position)<r &&
           i.triangleIndex != j.triangleIndex)) // ignore the case where i & j are the same point
           return vec3(0, 0, 0); // return black
    }
    return (lightColor * max(glm::dot(r_hat, n_hat), 0.0f)) / (4.0f * PI * r*r);
}

void emitPhotons(int nPhotons)
{
    float x,y,z;
    vec3 photon_normal;

    for(int i=0; i<nPhotons; i++)
    {
        do { // use simple rejection sampling to find diffuse photon direction
            x = random_num();
            y = random_num();
            z = random_num();
        } while ( x*x + y*y + z*z > 1 );

        photon_normal = vec3(x,y,z);
        Photon p(photon_normal, lightPos, lightPower/(float) nPhotons);

        trace_photon(p);
    }
}

void trace_photon(Photon& p)
{
    Intersection i, j;
    vec3 x0, x1; float t0,t1;
   
    if(!ClosestIntersection(p.getSource(), p.getNormal(), triangles, spheres, i)) return;
    
    if(i.sphereIndex >= 0)
    {
        // j is the intersection after refraction
        refract(spheres[i.sphereIndex], p.getNormal(), triangles, spheres, i, j);

        p.setDestination(j.position);
        photons.push_back(p);
        std::swap(i,j);
    }
    else
    {
        p.setDestination(i.position);
        if(p.getBounces()!=0)
            photons.push_back(p);
    }
   
    // New photon
    Photon p2(random_direction(triangles[i.triangleIndex].normal),   // normal vector
              i.position,                                                 // position of source
              p.getEnergy()*triangles[i.triangleIndex].color/(float) sqrt(p.getBounces()+1), // power
              p.getBounces()+1);                                          // increase bounces

    if (random_num(0) < 0.8)
        trace_photon(p2);
}

ivec2 VertexShader(vec3 p)
{
    ivec2 p2;
    vec3 p_prime = p-cameraPos;
    p2.x = focalLength*(p_prime.x/p_prime.z) + SCREEN_WIDTH/2;
    p2.y = focalLength*(p_prime.y/p_prime.z) + SCREEN_HEIGHT/2;

    return p2;
}

void DrawPhoton(Photon p)
{
    ivec2 p2 = VertexShader(p.getDestination());
    PutPixelSDL(screen, p2.x, p2.y, vec3(1,1,1));
}

void DrawPhotonPath(Photon p)
{
    ivec2 a = VertexShader(p.getSource());
    ivec2 b = VertexShader(p.getDestination());

    ivec2 delta = glm::abs(a - b);
    int pixels = glm::max(delta.x, delta.y) + 1;
    vector<ivec2> line(pixels);
    Interpolate(a, b, line);

    for(int l=0; l<line.size(); l++)
    {
        PutPixelSDL(screen, line[l].x, line[l].y, p.getEnergy());
    }
}

