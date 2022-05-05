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
vec3 lightColor = 2.f * vec3(1, 1, 1);
vec3 indirectLight = 0.1f*vec3( 1, 1, 1 );

int k = 600; // nearest photons TODO
float maxDist2 = 1.0; // TODO

// Intersection
struct Intersection
{
    vec3 position;
    float distance;
    int triangleIndex;
};

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

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    LoadTestModel(triangles);
    
    int nPhotons = 200000;
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

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    SDL_FillRect(screen, 0, 0);
    vec3 color, dir;
    Intersection closestIntersection;
	for(int y=0; y<SCREEN_HEIGHT; ++y)
	{
		for(int x=0; x<SCREEN_WIDTH; ++x)
		{
            vec3 dir(x-SCREEN_WIDTH/2, y-SCREEN_HEIGHT/2, focalLength);
            if(ClosestIntersection(R*cameraPos, dir, triangles,
                                   closestIntersection))
            {
                vector<int> neighbor_photons_idx = 
                              photonmap.searchKNearest(closestIntersection.position, 
                                                       k, maxDist2);

                color = vec3(0,0,0);
                for(int i=0; i<neighbor_photons_idx.size(); i++)
                    color += photons[neighbor_photons_idx[i]].getEnergy();
                color /= (float) k;

                PutPixelSDL(screen, x, y, triangles[closestIntersection.triangleIndex].color *
                            color+indirectLight);
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
        Photon p(photon_normal, lightPos, lightColor);

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
    photons.push_back(p);
    DrawPhoton(p);

    Photon p2(uniform_random_normal(triangles[i.triangleIndex].normal),
              i.position, triangles[i.triangleIndex].color/(float)(p.getBounces()+1),
              p.getBounces()+1);

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
    glm::vec2 step = glm::vec2(b-a) / float(max(N-1,1));
    glm::vec2 current( a );
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

    //printf("a(%d, %d)", a.x, a.y);
    //printf("b(%d, %d)\n", b.x, b.y);

    ivec2 delta = glm::abs(a - b);
    int pixels = glm::max(delta.x, delta.y) + 1;
    vector<ivec2> line(pixels);
    Interpolate(a, b, line);

    for(int l=0; l<line.size(); l++)
    {
        PutPixelSDL(screen, line[l].x, line[l].y, p.getEnergy());
    }
}

