/**
 * @file photon_map.cpp
 * @author Georgios Hadjiantonis
 * @brief Photon-mapping and ray tracing.
 *
 */
#include "photon_map.h"

#include <iostream>

#include <glm/glm.hpp>

#include "shapes.h"
#include "test_model.h"

using namespace std;

using glm::ivec2;
using glm::mat3;
using glm::vec2;
using glm::vec3;

int main(int argc, char *argv[])
{
    srand(0); // Set random seed

    int t2;
    float dt;
    screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
    t = SDL_GetTicks(); // Set start value for timer.
    LoadTestModel(triangles, spheres);

    // Emit photons from light source and build the photon-map
    EmitPhotons(NUM_PHOTONS);
    photonmap.SetPhotons(photons.data(), photons.size());
    photonmap.BuildTree();

    // Compute photon-map creation time
    t2 = SDL_GetTicks();
    dt = float(t2 - t);
    t = t2;
    cout << "Photon-map creation time: " << dt << " ms." << endl;
    cout << "Photon-map size: " << photons.size() << '\n';

    // Render scene
    RayTrace();

    // Compute rendering time
    t2 = SDL_GetTicks();
    dt = float(t2 - t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;

    // Save rendered image
    SDL_SaveBMP(screen, "screenshot.bmp");

    return 0;
}

vec3 GetRadianceTriangle(const Intersection &i)
{
    vec3 color = vec3(0, 0, 0);
    vec3 delta_phi;
    float wpc;         // weight
    float dp;          // distance between photon and intersection point
    float r_sqr = 0.0; // max sqr distance from k-th photon

    // Get k-nearest photons
    vector<NeighborPhoton> neighbor_photons = photonmap.SearchKNearest(i.position, K_NEAREST, r_sqr);
    for (int p = 0; p < neighbor_photons.size(); p++)
    {
        // Cone-filter
        dp = neighbor_photons[p].dist;
        wpc = 1 - dp / (CONE_FILTER_CONST * sqrt(r_sqr));

        // Photon power
        delta_phi = max(glm::dot(-photons[neighbor_photons[p].index].direction,
                                 triangles[i.triangle_index].normal),
                        0.0f) *
                    photons[neighbor_photons[p].index].energy;
        color += wpc * delta_phi;
    }
    color /= (1 - 2 / (3 * CONE_FILTER_CONST) * PI * r_sqr);

    color += DirectLight(i); // Direct light
    color *= triangles[i.triangle_index].color;

    return color;
}

vec3 GetRadianceSphere(const Intersection &i)
{
    Intersection intersection_refract, intersection_reflect;

    Refract(spheres[i.sphere_index], glm::normalize(i.position - CAMERA_POS),
            triangles, spheres, i, intersection_refract);

    float f = Fresnel(glm::normalize(i.position - CAMERA_POS),
                      glm::normalize(i.position - spheres[i.sphere_index].center));
    float c = max(0.0, min(1.0, 1.5 * f));

    vec3 dir_reflect = Reflect(glm::normalize(i.position - CAMERA_POS),
                               glm::normalize(i.position - spheres[i.sphere_index].center));
    ClosestIntersection(i.position, dir_reflect, triangles, spheres, intersection_reflect);

    return (1 - c) * GetRadianceTriangle(intersection_refract) +
           (c)*GetRadianceTriangle(intersection_reflect);
}

void RayTrace()
{
    if (SDL_MUSTLOCK(screen))
        SDL_LockSurface(screen);

    SDL_FillRect(screen, 0, 0);
    Intersection closest_intersection;

    for (int y = 0; y < SCREEN_HEIGHT; ++y)
        for (int x = 0; x < SCREEN_WIDTH; ++x)
        {
            vec3 dir(x - SCREEN_WIDTH / 2, y - SCREEN_HEIGHT / 2, FOCAL_LENGTH);

            if (!ClosestIntersection(CAMERA_POS, dir, triangles, spheres, closest_intersection))
                continue;

            if (closest_intersection.sphere_index >= 0)
            {
                PutPixelSDL(screen, x, y, GetRadianceSphere(closest_intersection));
            }
            else
            {
                PutPixelSDL(screen, x, y, GetRadianceTriangle(closest_intersection));
            }
        }

    if (SDL_MUSTLOCK(screen))
        SDL_UnlockSurface(screen);

    SDL_UpdateRect(screen, 0, 0, 0, 0);
}

vec3 DirectLight(const Intersection &i)
{
    // Distance of object from light source
    float r = glm::distance(LIGHT_POS, i.position);

    // Direction from surface to the light source
    vec3 r_hat = glm::normalize(LIGHT_POS - i.position);

    // Normal pointing out from the surface
    vec3 n_hat = triangles[i.triangle_index].normal;

    Intersection j;
    if (ClosestIntersection(LIGHT_POS, glm::normalize(i.position - LIGHT_POS),
                            triangles, spheres, j))
    {
        if (j.sphere_index >= 0 ||
            (glm::distance(LIGHT_POS, j.position) < r &&
             i.triangle_index != j.triangle_index)) // ignore the case where i & j are the same point
            return vec3(0, 0, 0);                   // return black
    }
    return (LIGHT_COLOR * max(glm::dot(r_hat, n_hat), 0.0f)) / (4.0f * PI * r * r);
}

void EmitPhotons(int num_photons)
{
    float x, y, z;
    vec3 photon_dir;

    for (int i = 0; i < num_photons; i++)
    {
        do
        { // Use rejection sampling to find photon direction
            x = GetRandomNum();
            y = GetRandomNum();
            z = GetRandomNum();
        } while (x * x + y * y + z * z > 1);

        photon_dir = vec3(x, y, z);

        Photon p(photon_dir, LIGHT_POS, LIGHT_POWER / (float)num_photons);
        TracePhoton(p);
    }
}

void TracePhoton(Photon &p)
{
    Intersection i, j;
    vec3 x0, x1;
    float t0, t1;

    if (!ClosestIntersection(p.source, p.direction, triangles, spheres, i))
        return;

    if (i.sphere_index >= 0)
    {
        // j is the intersection after refraction
        Refract(spheres[i.sphere_index], p.direction, triangles, spheres, i, j);

        p.destination = j.position;
        photons.push_back(p);
        std::swap(i, j);
    }
    else
    {
        p.destination = i.position;
        if (p.bounces != 0)
            photons.push_back(p);
    }

    // New photon
    Photon p2(GetRandomDirection(triangles[i.triangle_index].normal),                    // normal vector
              i.position,                                                                // position of source
              p.energy * triangles[i.triangle_index].color / (float)sqrt(p.bounces + 1), // power
              p.bounces + 1);                                                            // increase bounces

    if (GetRandomNum(0) < 0.8)
        TracePhoton(p2);
}

ivec2 VertexShader(vec3 p)
{
    ivec2 p2;
    vec3 p_prime = p - CAMERA_POS;
    p2.x = FOCAL_LENGTH * (p_prime.x / p_prime.z) + SCREEN_WIDTH / 2;
    p2.y = FOCAL_LENGTH * (p_prime.y / p_prime.z) + SCREEN_HEIGHT / 2;

    return p2;
}

void DrawPhoton(Photon p)
{
    ivec2 p2 = VertexShader(p.destination);
    PutPixelSDL(screen, p2.x, p2.y, vec3(1, 1, 1));
}

void DrawPhotonPath(Photon p)
{
    ivec2 a = VertexShader(p.source);
    ivec2 b = VertexShader(p.destination);

    ivec2 delta = glm::abs(a - b);
    int pixels = glm::max(delta.x, delta.y) + 1;
    vector<ivec2> line(pixels);
    Interpolate(a, b, line);

    for (int l = 0; l < line.size(); l++)
        PutPixelSDL(screen, line[l].x, line[l].y, p.energy);
}
