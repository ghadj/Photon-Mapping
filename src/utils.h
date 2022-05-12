#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <random>
#include <glm/glm.hpp>

#include "shapes.h"

using namespace std;

using glm::ivec2;
using glm::vec2;
using glm::vec3;

const float PI = 3.14159265358979323846;

// Intersection
struct Intersection
{
    vec3 position;
    float distance;
    int sphere_index;
    int triangle_index;
};

// Backface culling
bool IsVisible(vec3 n, vec3 s)
{
    return glm::dot(n, s) <= 0.0;
}

bool ClosestIntersection(const vec3 start, const vec3 dir,
                         const vector<Triangle> &triangles,
                         const vector<Sphere> &spheres,
                         Intersection &intersection)
{
    intersection.distance = std::numeric_limits<float>::max();
    intersection.triangle_index = -1;
    intersection.sphere_index = -1;

    bool intersects = false;
    float t0, t1; // distance
    vec3 x0, x1;  // intersection points

    for (int t_i = 0; t_i < triangles.size(); t_i++)
    {
        // Backface culling
        if (!IsVisible(triangles[t_i].normal, dir))
        {
            continue;
        }
        if (triangles[t_i].Intersect(start, dir, x0, t0) &&
            t0 < intersection.distance)
        { // valid intersection
            intersects = true;
            intersection.position = x0;
            intersection.distance = t0;
            intersection.triangle_index = t_i;
        }
    }

    // Check intersection with spheres
    for (int s_i = 0; s_i < spheres.size(); s_i++)
    {
        if (spheres[s_i].Intersect(start, dir, x0, x1, t0, t1) &&
            ((t0 > 0.0001 && t0 < intersection.distance) ||
             (t1 > 0.0001 && t1 < intersection.distance)))
        {
            intersects = true;
            intersection.position = t0 > 0 ? x0 : x1;
            intersection.distance = t0 > 0 ? t0 : t1;
            intersection.sphere_index = s_i;
        }
    }
    return intersects;
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2> &result)
{
    int N = result.size();
    vec2 step = vec2(b - a) / float(max(N - 1, 1));
    vec2 current(a);

    for (int i = 0; i < N; i++)
    {
        result[i] = current;
        current += step;
    }
}

float GetRandomNum(float min = -1.0, float max = 1.0)
{
    float r = (float)rand() / (float)RAND_MAX;
    return min + r * (max - min);
}

vec3 GetRandomDirection(const vec3 &n)
{
    float z = sqrt(GetRandomNum(0));
    float r = sqrt(1.0 - z * z);
    float phi = 2.0 * PI * GetRandomNum(0);
    float x = r * cos(phi);
    float y = r * sin(phi);

    // local orthogonal coordinate system around n
    vec3 w, u, v;
    w = n;
    u = glm::normalize(glm::cross(abs(w.x) > .1 ? vec3(0, 1, 0) : vec3(1, 0, 0), w));
    v = glm::cross(w, u);

    return x * u + y * v + z * w;
}

#endif
