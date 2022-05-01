#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const float PI = 3.14159265358979323846; // pi
const int SCREEN_WIDTH = 300;
const int SCREEN_HEIGHT = 300;
const float TRANSL_STEP = 0.1;
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
vec3 lightColor = 14.f * vec3(1, 1, 1);
vec3 indirectLight = 0.5f*vec3( 1, 1, 1 );

// Intersection
struct Intersection
{
    vec3 position;
    float distance;
    int triangleIndex;
};

// ----------------------------------------------------------------------------
// FUNCTIONS

bool ClosestIntersection(vec3, vec3, const vector<Triangle>&, Intersection&);
vec3 DirectLight(const Intersection&);
void Update();
void Draw();

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    LoadTestModel(triangles);
	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	//SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2-t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;
    Uint8* keystate = SDL_GetKeyState( 0 );

    // 4. Moving the camera
    if( keystate[SDLK_UP] )
    {
        cameraPos[2] += TRANSL_STEP;
    }
    if( keystate[SDLK_DOWN] )
    {
        cameraPos[2] -= TRANSL_STEP;
    }
    if( keystate[SDLK_LEFT] )
    {
        yaw -= TRANSL_STEP;
    }
    if( keystate[SDLK_RIGHT] )
    {
        yaw += TRANSL_STEP;
    }

    // Set camera rotation matrix
    R = mat3(cos(yaw),  0, sin(yaw),
                0,      1,    0,   
             -sin(yaw), 0, cos(yaw));

    // 5. Illumination
    if( keystate[SDLK_w] )
    {
        lightPos.y -= TRANSL_STEP;
    }
    if( keystate[SDLK_s] )
    {
        lightPos.y += TRANSL_STEP;
    }
    if( keystate[SDLK_d] )
    {
        lightPos.x += TRANSL_STEP;
    }
    if( keystate[SDLK_a] )
    {
        lightPos.x -= TRANSL_STEP;
    }
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
                // 3. Tracing Rays
                //color = triangles[closestIntersection.triangleIndex].color;
			    
                // 5. Illumination (Figure 4)
                //color = DirectLight(closestIntersection);

                // 5. Illumination (Figure 6)
                //color = triangles[closestIntersection.triangleIndex].color * 
                //        DirectLight(closestIntersection);

                // 5.2 Indirect Illumation (Figure 8)
                color = triangles[closestIntersection.triangleIndex].color * 
                        (DirectLight(closestIntersection) + indirectLight);

                PutPixelSDL(screen, x, y, color);
            }
		}
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle>& triangles,
                         Intersection& closestIntersection)
{
    closestIntersection.distance = std::numeric_limits<float>::max();
    bool intersect = false;
    vec3 v0, v1,v2, e1, e2, b, x;

    for(int t_i=0; t_i<triangles.size(); t_i++)
    {
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
