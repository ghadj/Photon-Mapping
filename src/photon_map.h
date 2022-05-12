/**
 * @file photon_map.h
 * @author Georgios Hadjiantonis
 * @brief Photon-mapping and ray tracing.
 *
 */
#ifndef PHOTON_MAP_H
#define PHOTON_MAP_H

#include <glm/glm.hpp>
#include <SDL.h>
#include "libs/SDLauxiliary.h"

#include "kdtree.h"
#include "photon_tracer.h"
#include "utils.h"

using glm::ivec2;
using glm::vec3;

// Screen Dimensions
const int SCREEN_WIDTH = 300;
const int SCREEN_HEIGHT = 300;

// Camera parameters
const float FOCAL_LENGTH = SCREEN_HEIGHT;
const vec3 CAMERA_POS(0, 0, -3); // Position of the camera

// Light parameters
const vec3 LIGHT_POS(0, -0.5, -0.4);
const vec3 LIGHT_COLOR = 3.0f * vec3(1, 1, 0.95);
const vec3 LIGHT_POWER = 60.f * vec3(1, 1, 1);

// Photon-map parameters
const int K_NEAREST = 500;           // Number of nearest photons
const float CONE_FILTER_CONST = 1.2; // Cone-filter constant
const int NUM_PHOTONS = 5000;        // Total number of photons from light source

// Global variables
vector<Photon> photons;
KdTree photonmap;

vector<Triangle> triangles; // All the triangles of the scene
vector<Sphere> spheres;     // All the spheres of the scene

SDL_Surface *screen;
int t; // Time

/**
 * @brief
 *
 * @return vec3
 */
vec3 DirectLight(const Intersection &);

/**
 * @brief Renders the scene by ray-tracing and utilizing the created photon-map.
 *
 */
void RayTrace();

/**
 * @brief Emission of photons from a diffuse point light.
 *
 * @param num_photons Number of photons to be emitted.
 */
void EmitPhotons(int num_photons);

/**
 * @brief Traces the given photon in the scene.
 *
 * Photons are added at the end to the global photon vector on their second
 * bounce. (Direct lighting is calculated separately) Only exception is for the
 * caustics, when the photon is refracted/reflected by a (glass) sphere. 
 *
 * The power of a photon is decreased by a factor of the sqrt of the number of
 * bounces and color of the material. 
 *
 * For simplicity all triangles are considered as ideally diffused surfaces, 
 * i.e., the photon has equal probability to be reflected on all directions.
 *
 * Given a predefined probability the photon is either reflected or terminated.
 *
 * @param p A photon.
 */
void TracePhoton(Photon &p);

/**
 * @brief Draws the position of the given photon on the screen.
 *
 * @param p A photon.
 */
void DrawPhoton(Photon p);

/**
 * @brief Draws the path of the given photon on the screen, by interpolating
 *        origin and destination coordinates.
 *
 * @param p A photon.
 */
void DrawPhotonPath(Photon p);

/**
 * @brief Transforms 3D coordinates to 2D given the parameters of the camera.
 *
 * @param p      3D coordinates.
 * @return ivec2 Tranformed 2D coordinates.
 */
ivec2 VertexShader(vec3 p);

/**
 * @brief Calculates the radiance at the given intersection point on a triangle.
 *
 * Implementation of eq. (11) using a cone filter from: 
 * Jensen, H. W., & Christensen, N. J. (2000). A practical guide to global 
 * illumination using photon maps. SIGGRAPH 2000
 *
 * @param i     Intersection point on a triangle.
 * @return vec3 Radiance/color in RGB form.
 */
vec3 GetRadianceTriangle(const Intersection &i);

/**
 * @brief Calculates the radiance at the given intersection point on a sphere.
 *
 * Used an approximation of the Fresnel coefficient considering both reflected
 * and refracted light.
 *
 * @param i     Intersection point on a sphere.
 * @return vec3 Radiance/color in RGB form.
 */
vec3 GetRadianceSphere(const Intersection &i);

#endif
