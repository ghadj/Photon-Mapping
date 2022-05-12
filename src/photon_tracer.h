#ifndef PHOTON_TRACER_H
#define PHOTON_TRACER_H

#include <vector>

#include <glm/glm.hpp>

#include "shapes.h"
#include "utils.h"

using glm::vec3;

struct Photon
{
    vec3 direction;
    vec3 source;
    vec3 destination;
    vec3 energy; // incoming photon power
    int bounces;

    Photon(vec3 direction, vec3 source, vec3 energy, int bounces = 0)
        : direction(direction), source(source), energy(energy), bounces(bounces) {}
};

// Reference. Based on
// https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/
// reflection-refraction-fresnel
bool Refract(const vec3 &photon_dir,     // direction of photon
             const vec3 &surface_normal, // surface normal
             vec3 &T,                    // refracted direction
             const float &ior = 1.75)    // index of refraction
{
    vec3 N = surface_normal;
    float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(photon_dir, N));
    float n1 = 1;   // refractive index being left
    float n2 = ior; // refractive index being entered
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

    n = (n1 / n2);
    c = 1 - n * n * (1 - cosi * cosi);

    if (c < 0)
        return false;

    T = n * photon_dir + (n * cosi - sqrtf(c)) * N;
    return true;
}

bool Refract(const Sphere s, const vec3 photon_dir,
             const vector<Triangle> &triangles,
             const vector<Sphere> &spheres,
             const Intersection &i, // Intersection point with sphere
             Intersection &j) // Intersection point after refraction
{
    vec3 sphere_normal = glm::normalize(i.position - s.center); // normal at intersection point
    vec3 Ti, Tj;                                                // refracted direction

    // Refract from outside
    if (!Refract(photon_dir, sphere_normal, Ti))
        return false;

    // Find intersection from inside
    vec3 x0, x1;
    float t0, t1;
    s.Intersect(i.position, Ti, x0, x1, t0, t1);
    sphere_normal = glm::normalize(x1 - s.center);

    // Refract to ouside
    if (!Refract(Ti, sphere_normal, Tj))
        return false;
    // Find intersection ouside
    ClosestIntersection(x1, Tj, triangles, spheres, j);
    return true;
}

vec3 Reflect(vec3 photon_dir, vec3 surface_normal)
{
    return photon_dir - 2.0f * surface_normal * glm::dot(photon_dir, surface_normal);
}

// Reference.
// https://developer.download.nvidia.com/SDK/9.5/Samples/DEMOS/Direct3D9/src/
// HLSL_FresnelReflection/docs/FresnelReflection.pdf
//
// Schlick’s Fresnel approximation function
float Fresnel(vec3 photon_dir, vec3 surface_normal, const float &ior = 1.75)
{ // light and normal are assumed to be normalized
    float r0 = pow(1.0 - ior, 2.0) / pow(1.0 + ior, 2.0);
    float cosx = -glm::dot(surface_normal, photon_dir);

    return r0 + (1.0 - r0) * pow(1.0 - cosx, 5.0);
}

#endif
