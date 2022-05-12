/**
 * @file   photon_tracer.h
 * @author Georgios Hadjiantonis
 * @brief  Description of photon instance and related functions for its
 *         reflection & refraction. 
 */
#ifndef PHOTON_TRACER_H
#define PHOTON_TRACER_H

#include <vector>

#include <glm/glm.hpp>

#include "shapes.h"
#include "utils.h"

using glm::vec3;

/**
 * @brief Photon instance.
 *
 */
struct Photon
{
    vec3 direction;   // direction from source
    vec3 source;      // coordinates from where the photon originated
    vec3 destination; // where the photon intersected with a surface
    vec3 energy;      // incoming photon power
    int bounces;      // number of bounces at the current state

    Photon(vec3 direction, vec3 source, vec3 energy, int bounces = 0)
        : direction(direction), source(source), energy(energy), bounces(bounces) {}
};

/**
 * @brief Calculates the direction of the photon after refraction with a
 *        surface.
 *
 * Reference:
 * https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/
 * reflection-refraction-fresnel
 *
 * @param photon_dir     Direction of the photon.
 * @param surface_normal Normal vector of the surface.
 * @param T              Direction of the photon after refraction.
 * @param ior            Index of refraction (depends on the material).
 * @return true          The photon hits the surface and is refracted. 
 * @return false         There is no refraction in this case (total internal reflection). 
 */
bool Refract(const vec3 &photon_dir, const vec3 &surface_normal, vec3 &T,
             const float &ior = 1.75)
{
    vec3 N = surface_normal;
    float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(photon_dir, N));
    float n1 = 1;   // refractive index of origin
    float n2 = ior; // refractive index of destination
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

/**
 * @brief Calculates the intersection of a photon after hitting and being
 *        refracted with the given sphere.
 *
 * @param s          The sphere the photon hits.
 * @param photon_dir Direction of the photon.
 * @param triangles  Vector of all the triangles in the scene.
 * @param spheres    Vector of all the spheres in the scene.
 * @param i          Intersection with the given sphere.
 * @param j          Intersection after the photon is refracted by the sphere.
 * @return true      The photon hits the sphere, is refracted and hits another
 *                   surface succesfully.
 * @return false     There is no refraction (total internal reflection) or the
 *                   photon does not intersect with a surface after refraction.
 */
bool Refract(const Sphere s, const vec3 photon_dir,
             const vector<Triangle> &triangles, const vector<Sphere> &spheres,
             const Intersection &i, Intersection &j) 
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

/**
 * @brief Calculates the reflection vector.
 *
 * @param photon_dir     Direction of the photon.
 * @param surface_normal Normal vector of the surface.
 * @return vec3          Direction after reflection.
 */
vec3 Reflect(vec3 photon_dir, vec3 surface_normal)
{
    return photon_dir - 2.0f * surface_normal * glm::dot(photon_dir, surface_normal);
}

/**
 * @brief Implementation of the Schlick’s Fresnel approximation function.
 *
 * Assumes that direction vector and surface normal to be normalized.
 *
 * Reference:
 * https://developer.download.nvidia.com/SDK/9.5/Samples/DEMOS/Direct3D9/src/
 * HLSL_FresnelReflection/docs/FresnelReflection.pdf
 *
 * @param photon_dir     Direction of the photon.
 * @param surface_normal Normal vector of the surface.
 * @param ior            Index of refraction (depends on the material).
 * @return float         Fresnel coefficient value.
 */
float Fresnel(vec3 photon_dir, vec3 surface_normal, const float &ior = 1.75)
{ 
    float r0 = pow(1.0 - ior, 2.0) / pow(1.0 + ior, 2.0);
    float cosx = -glm::dot(surface_normal, photon_dir);

    return r0 + (1.0 - r0) * pow(1.0 - cosx, 5.0);
}

#endif
