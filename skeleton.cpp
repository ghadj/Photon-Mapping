#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <random>
#include "SDLauxiliary.h"
#include "TestModel.h"

#include "photon.h"
#include "kdtree.h"

using namespace std;

using glm::vec3;
using glm::mat3;
using glm::ivec2;
using glm::vec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const float PI = 3.14159265358979323846; // pi
const int SCREEN_WIDTH = 300;
const int SCREEN_HEIGHT = 300;
SDL_Surface* screen;
vector<Triangle> triangles; // all the trianlges of the scene
int t; // time

// Camera
float focalLength = SCREEN_HEIGHT;
vec3 cameraPos(0, 0, -3); // position of the camera
mat3 R; // rotation matrix of the camera
float yaw; // angle of the camera around y axis 

// Light
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColor = 7.f * vec3(1, 1, 1);
vec3 lightPower = 500.f * vec3(1, 1, 1);
vec3 indirectLight = 0.1f*vec3( 1, 1, 1 );

int k = 700; // nearest photons TODO
float filter_const = 1;

// Intersection
struct Intersection
{
    vec3 position;
    float distance;
    int triangleIndex;
};

int nPhotons = 200000;
vector<Photon> photons;
KdTree photonmap;

// ----------------------------------------------------------------------------
// FUNCTIONS

bool ClosestIntersection(vec3, vec3, const vector<Triangle>&, Intersection&);
vec3 DirectLight(const Intersection&);
void Draw();
void emitPhotons(int nPhotons);
void trace_photon(Photon& p);
void DrawPhoton(Photon p);
void DrawPhotonPath(Photon p);
float random_num(float min=-1, float max=1);
bool isVisible(vec3 n, vec3 s);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    LoadTestModel(triangles);
    
    emitPhotons(nPhotons);
    photonmap.setPhotons(photons.data(), photons.size());
    photonmap.buildTree();
    cout<<"Photon-map size: "<<photons.size()<<'\n';

    Draw();
    
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
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
    vector<NeighborPhoton> neighbor_photons = photonmap.searchKNearest(i.position, k, r_sqr);

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

    // TODO remove
    //cout<<"color: ("<<color.x<<','<<color.y<<','<<color.z<<")\n";

    return color;
}

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    SDL_FillRect(screen, 0, 0);
    Intersection closestIntersection;
	for(int y=0; y<SCREEN_HEIGHT; ++y)
	{
		for(int x=0; x<SCREEN_WIDTH; ++x)
		{
            vec3 dir(x-SCREEN_WIDTH/2, y-SCREEN_HEIGHT/2, focalLength);
            if(ClosestIntersection(R*cameraPos, dir, triangles,
                                   closestIntersection))
            {
                PutPixelSDL(screen, x, y, getRadianceEstimate(closestIntersection));
            }
		}
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

// Backface culling
bool isVisible(vec3 n, vec3 s)
{
    return glm::dot(n, s)<=0.0;
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles,
                         Intersection& closestIntersection)
{
    closestIntersection.distance = std::numeric_limits<float>::max();
    bool intersect = false;
    vec3 v0, v1,v2, e1, e2, b, x;

    for(int t_i=0; t_i<triangles.size(); t_i++)
    {
        if(!isVisible(triangles[t_i].normal, dir))
            continue;

        v0 = triangles[t_i].v0;
        v1 = triangles[t_i].v1;
        v2 = triangles[t_i].v2;
        e1 = v1 - v0;
        e2 = v2 - v0;
        b = start - v0;
        mat3 A( -dir, e1, e2 );
        x = glm::inverse( A ) * b; // x:t, y:u, z:v

        // Add equallity so that it includes points on the edges of the triang.
        if(x.x>=0 && x.y>=0 && x.z>=0 && (x.y+x.z)<=1 && 
           x.x<closestIntersection.distance)
        {   // valid intersection
            intersect = true;
            closestIntersection.position = v0 + x.y*e1 + x.z*e2;
            closestIntersection.distance = x.x;
            closestIntersection.triangleIndex = t_i;
        }
    }
    return intersect;
}

vec3 DirectLight(const Intersection& i)
{
    // Distance of object from light source
    float r = glm::distance(lightPos, i.position); 

    // Direction from surface to the light source
    // Reference: https://math.hws.edu/graphicsbook/c7/s2.html
    vec3 r_hat = glm::normalize(lightPos - i.position); 

    // Normal pointing out from the surface
    vec3 n_hat = triangles[i.triangleIndex].normal; 
    
    // 5. Illumation (Figures 4 & 6)
    //return (lightColor*max(glm::dot(r_hat, n_hat), 0.0f))/(4.0f*PI*r*r);

    // 5.1 Direct shadows (Figure 7 & 8)
    Intersection j;
    if(ClosestIntersection(lightPos, glm::normalize(i.position - lightPos), 
                           triangles, j))
    {
       if(glm::distance(lightPos, j.position)<r &&
          i.triangleIndex != j.triangleIndex) // ignore the case where i & j are 
                                              // the same point
           return vec3(0, 0, 0); // return black
    }
    return (lightColor*max(glm::dot(r_hat, n_hat), 0.0f))/(4.0f*PI*r*r);
}

float random_num(float min, float max)
{
	float r = (float)rand() / (float)RAND_MAX;
	return min + r * (max - min);
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

vec3 uniform_random_normal(const vec3& n) 
{
    float z = sqrt(random_num(0));
    float r = sqrt(1.0 - z * z);
    float phi = 2.0 * PI * random_num(0);
    float x = r * cos(phi);
    float y = r * sin(phi);

    // local orthogonal coordinate system around n
    vec3 w = n;
    vec3 u = glm::normalize(glm::cross(abs(w.x)>.1 ? vec3(0,1,0) : 
                            vec3(1,0,0), w));
    vec3 v = glm::cross(w, u);

    return x*u + y*v + z*w;
}

void trace_photon(Photon& p)
{
    Intersection i;
   
    if(!ClosestIntersection(p.getSource(), p.getNormal(), triangles, i))
    {
        return;
    }

    p.setDestination(i.position);

    if(p.getBounces()!=0)
        photons.push_back(p);
    
    //DrawPhoton(p);

    // New photon
    Photon p2(uniform_random_normal(triangles[i.triangleIndex].normal),   // normal vector
              i.position,                                                 // position of source
              p.getEnergy()*triangles[i.triangleIndex].color/(float) sqrt(p.getBounces()+1), // power
              p.getBounces()+1);                                          // increase bounces

    if (random_num(0)<0.8)
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
    PutPixelSDL(screen, p2.x, p2.y, p.getEnergy());
}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
    int N = result.size();
    vec2 step = vec2(b-a) / float(max(N-1,1));
    vec2 current( a );
    for( int i=0; i<N; ++i )
    {
        result[i] = current;
        current += step;
    }
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

