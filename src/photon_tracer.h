#ifndef PHOTON_TRACER_H
#define PHOTON_TRACER_H

#include <vector>
#include <glm/glm.hpp>

#include "shapes.h"
#include "utils.h"

using glm::vec3;

class Photon
{
private:
    vec3 normal;
    vec3 source;
    vec3 destination;
    vec3 energy; // incoming photon power
    int bounces;
     
public:
    Photon(vec3 normal, vec3 source, vec3 energy, int bounces=0) 
        : normal(normal), source(source), energy(energy), bounces(bounces){}

    vec3 getNormal() const
    {
        return this->normal;
    }

    vec3 getSource() const
    {
        return this->source;
    }

    void setDestination(vec3 d)
    {
        this->destination = d;
    }

    vec3 getDestination() const
    {
        return this->destination;
    }

    vec3 getEnergy() const
    {
        return this->energy;
    }
    
    int getBounces() const
    {
        return this->bounces;
    }

};

// Reference. Based on
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/
// reflection-refraction-fresnel
bool refract(const vec3& light_dir,      // direction of photon
             const vec3& surface_normal, // surface normal
             vec3& T,                    // refracted direction
             const float& ior=1.75)      // index of refraction
{ 
    vec3 N = surface_normal; 
    float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(light_dir, N)); 
    float n1 = 1;    // refractive index being left
    float n2 = ior;  // refractive index being entered
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

bool refract(Sphere s, vec3 photon_dir, 
             vector<Triangle>& triangles, 
             vector<Sphere> spheres,
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
    s.intersect(i.position, Ti, x0, x1, t0, t1);
    sphere_normal = glm::normalize(x1 - s.center);

    // Refract to ouside
    if(!refract(Ti, sphere_normal, Tj)) return false;
    // Find intersection ouside
    ClosestIntersection(x1, Tj, triangles, spheres, j);
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
float fresnel(vec3 light_dir, vec3 surface_normal, const float& ior=1.75)
{   // light and normal are assumed to be normalized
    float r0 = pow(1.0-ior, 2.0) / pow(1.0+ior, 2.0);  
    float cosx = -glm::dot(surface_normal, light_dir);

    return r0 + (1.0-r0) * pow(1.0 - cosx, 5.0);
} 

#endif
