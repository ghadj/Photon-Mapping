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
const int SCREEN_WIDTH = 100;
const int SCREEN_HEIGHT = 100;
SDL_Surface* screen;
vector<Triangle> triangles; // all the trianlges of the scene
vector<Sphere> spheres;
int t; // time

// Camera
float focalLength = SCREEN_HEIGHT;
vec3 cameraPos(0, 0, -3); // position of the camera
mat3 R; // rotation matrix of the camera
float yaw; // angle of the camera around y axis 

// Light
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColor = 4.0f * vec3(1, 1, 0.95);
vec3 lightPower = 400.f * vec3(1, 1, 1);

int k = 600; // nearest photons
float filter_const = 1;

// Intersection
struct Intersection
{
    vec3 position;
    float distance;
    int triangleIndex;
    int sphereIndex;
};

int nPhotons = 100000;
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
bool intersectSphere(const Sphere& s, const vec3& start, const vec3& dir,
                     vec3& x0, vec3& x1,   // intersection points
                     float& t0, float& t1); // distance from origin
ivec2 VertexShader(vec3 p);
bool refractPhoton(Sphere s, vec3 photon_dir, 
                   Intersection& i, // Intersection point with sphere
                   Intersection& j); // Intersection point after refraction
vec3 reflect(vec3 light_dir, vec3 surface_normal);
float fresnel(vec3 light_dir, vec3 surface_normal, const float& ior=1.66);

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

    Draw();
    
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

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    SDL_FillRect(screen, 0, 0);
    Intersection closestIntersection, intersection_refract, intersection_reflect;

	for(int y=0; y<SCREEN_HEIGHT; ++y)
		for(int x=0; x<SCREEN_WIDTH; ++x)
		{
            vec3 dir(x-SCREEN_WIDTH/2, y-SCREEN_HEIGHT/2, focalLength);
            if(ClosestIntersection(cameraPos, dir, triangles, closestIntersection))
            {
                if(closestIntersection.sphereIndex >= 0)
                {
                    refractPhoton(spheres[closestIntersection.sphereIndex],
                                  glm::normalize(closestIntersection.position-cameraPos), 
                                  closestIntersection, intersection_refract);
                    
                    float f = fresnel(glm::normalize(closestIntersection.position-cameraPos),
                                      glm::normalize(closestIntersection.position - 
                                                     spheres[closestIntersection.sphereIndex].center));

                    vec3 dir_re = reflect(glm::normalize(closestIntersection.position-cameraPos),
                                      glm::normalize(closestIntersection.position - 
                                                     spheres[closestIntersection.sphereIndex].center));
                    ClosestIntersection(closestIntersection.position, dir_re, triangles, intersection_reflect);

                    PutPixelSDL(screen, x, y, (1-f)*getRadianceEstimate(intersection_refract)+
                                            (f)*getRadianceEstimate(intersection_reflect));
                }
                else
                    PutPixelSDL(screen, x, y, getRadianceEstimate(closestIntersection));
            }
		}

	if(SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

// Backface culling
bool isVisible(vec3 n, vec3 s)
{
    return glm::dot(n, s) <= 0.0;
}

// Reference. Based on
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/
// reflection-refraction-fresnel
bool refract(const vec3& light_dir,      // direction of photon
             const vec3& surface_normal, // surface normal
             vec3& T,                    // refracted direction
             const float& ior=1.66)      // index of refraction
{ 
    vec3 N = surface_normal; 
    float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(light_dir, N)); 
    float n1 = 1;
    float n2 = ior; 
    float n, c;
    
    if (cosi < 0) 
    {
        cosi *= -1; 
    }
    else 
    {
        std::swap(n1, n2); 
        N *= -1; 
    }

    n = (n1/n2);
    c = 1 - n*n * (1 - cosi * cosi); 

    if(c<0) return false; 

    T = n * light_dir + (n * cosi - sqrtf(c)) * N; 
    return true;
} 

bool refractPhoton(Sphere s, vec3 photon_dir, 
                   Intersection& i, // Intersection point with sphere
                   Intersection& j) // Intersection point after refraction
{
    vec3 sphere_normal = glm::normalize(i.position - s.center); // normal at intersection point
    vec3 Ti, Tj; // refracted direction
    
    // Refract from outside
    if(!refract(photon_dir, sphere_normal, Ti)) return false;

    // Find intersection from inside
    vec3 x0, x1;
    float t0, t1;
    intersectSphere(s, i.position, Ti, x0, x1, t0, t1);
    sphere_normal = glm::normalize(x1 - s.center);

    // Refract to ouside
    if(!refract(Ti, sphere_normal, Tj)) return false;
    // Find intersection ouside
    ClosestIntersection(x1, Tj, triangles, j);
    return true;
}

vec3 reflect(vec3 light_dir, vec3 surface_normal)
{
    return light_dir - 2.0f * surface_normal * glm::dot(surface_normal, light_dir);
}

// Reference.
// https://developer.download.nvidia.com/SDK/9.5/Samples/DEMOS/Direct3D9/src/
// HLSL_FresnelReflection/docs/FresnelReflection.pdf
//
// Schlick’s Fresnel approximation function
float fresnel(vec3 light_dir, vec3 surface_normal, const float& ior)
{   // light and normal are assumed to be normalized
    float r0 = pow(1.0-ior, 2.0) / pow(1.0+ior, 2.0);  
    float cosx = -glm::dot(surface_normal, light_dir);

    return r0 + (1.0-r0) * pow(1.0 - cosx, 5.0);
} 

bool intersectSphere(const Sphere& s, const vec3& start, const vec3& dir,
                     vec3& x0, vec3& x1,   // intersection points
                     float& t0, float& t1) // distance from origin
{ 
    vec3 center_trn = start- s.center; // translated center of sphere
    float A = glm::dot(dir, dir);   
    float B = 2.0 * glm::dot(dir, center_trn);   
    float C = glm::dot(center_trn, center_trn) - pow(s.radius, 2.0); 
    float Delta = pow(B, 2.0) - 4*A*C; 

    if(Delta == 0.0)      // there is only one solution
    {
        t0 = t1 = - 0.5 * B / A;
        x0 = start + t0*dir;
        x1 = start + t1*dir;
        return true;
    }
    else if (Delta > 0.0) // there are two solutions
    { 
        t0 = (-B + sqrt(Delta))/(2*A);
        t1 = (-B - sqrt(Delta))/(2*A);

        if(t0<0 && t1<0) return false;

        // Note: now there is the possibility one of t0 or t1 to be negative.
        // In that case the origin of the ray is inside the sphere.
    
        // Swap so that t1 corresponds to the farthest intersection point
        if (t0 > t1) std::swap(t0, t1); 

        x0 = start + t0*dir;
        x1 = start + t1*dir;

        return true;
    }
    //cout<<"RETURN FALSE";
    return false;         // no intersection
}

bool intersectTriangle(vec3 v0, vec3 v1, vec3 v2, vec3 start, vec3 dir, 
                       vec3& x0, float& t0)
{
    vec3 e1, e2, b, x;
    e1 = v1 - v0;
    e2 = v2 - v0;
    b = start - v0;
    mat3 A(-dir, e1, e2);
    x = glm::inverse(A) * b; // x:t, y:u, z:v

    // Add equallity so that it includes points on the edges of the triang.
    if(x.x>=0 && x.y>=0 && x.z>=0 && (x.y+x.z)<=1)
    {
        x0 = v0 + x.y*e1 + x.z*e2;
        t0 = x.x;
        return true;
    }
    return false;
}

bool ClosestIntersection(const vec3 start, const vec3 dir,
                         const vector<Triangle>& triangles,
                         Intersection& closestIntersection)
{
    closestIntersection.distance = std::numeric_limits<float>::max();
    closestIntersection.triangleIndex = -1;
    closestIntersection.sphereIndex = -1;

    bool intersect = false;
    float t0, t1; // distance
    vec3 x0, x1;  // intersection points

    for(int t_i=0; t_i<triangles.size(); t_i++)
    {
        // Backface culling
        if(!isVisible(triangles[t_i].normal, dir))
            continue;

        if(intersectTriangle(triangles[t_i].v0, triangles[t_i].v1, triangles[t_i].v2,
                             start, dir, x0, t0) && t0<closestIntersection.distance)
        {   // valid intersection
            intersect = true;
            closestIntersection.position = x0;
            closestIntersection.distance = t0;
            closestIntersection.triangleIndex = t_i;
        }
    }

    // Check intersection with spheres
    for(int s_i=0; s_i<spheres.size(); s_i++)
    {
        if(intersectSphere(spheres[s_i], start, dir, x0, x1, t0, t1) &&
           ((t0>0.001 && t0<closestIntersection.distance) || 
            (t1>0.001 && t1<closestIntersection.distance) ))
        {
            intersect = true;
            closestIntersection.position = t0>0? x0 : x1;
            closestIntersection.distance = t0>0? t0 : t1;
            closestIntersection.sphereIndex = s_i;
        }
    }

    return intersect;
}

vec3 DirectLight(const Intersection& i)
{
    // Distance of object from light source
    float r = glm::distance(lightPos, i.position); 

    // Direction from surface to the light source
    vec3 r_hat = glm::normalize(lightPos - i.position); 

    // Normal pointing out from the surface
    vec3 n_hat = triangles[i.triangleIndex].normal; 
    
    // 5.1 Direct shadows (Figure 7 & 8)
    Intersection j;
    if(ClosestIntersection(lightPos, glm::normalize(i.position - lightPos), 
                           triangles, j))
    {
        if(j.sphereIndex >= 0 ||
           (glm::distance(lightPos, j.position)<r &&
           i.triangleIndex != j.triangleIndex)) // ignore the case where i & j are the same point
           return vec3(0, 0, 0); // return black
    }
    return (lightColor * max(glm::dot(r_hat, n_hat), 0.0f)) / (4.0f * PI * r*r);
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
    Intersection i, j;
    vec3 x0, x1; float t0,t1;
   
    if(!ClosestIntersection(p.getSource(), p.getNormal(), triangles, i)) return;
    
    if(i.sphereIndex >= 0)
    {
        // j is the intersection after refraction
        refractPhoton(spheres[i.sphereIndex], p.getNormal(), i, j);

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
    Photon p2(uniform_random_normal(triangles[i.triangleIndex].normal),   // normal vector
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

