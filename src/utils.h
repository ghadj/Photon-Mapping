#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <random>
#include <glm/glm.hpp>

#include "shapes.h"

using namespace std;

const float PI = 3.14159265358979323846;

// Intersection
struct Intersection
{
    glm::vec3 position;
    float distance;
    int sphereIndex;
    int triangleIndex;
};

// Backface culling
bool isVisible(glm::vec3 n, glm::vec3 s)
{
    return glm::dot(n, s) <= 0.0;
}

bool ClosestIntersection(const glm::vec3 start, const glm::vec3 dir,
                         const vector<Triangle>& triangles,
                         const vector<Sphere>& spheres,
                         Intersection& intersection)
{
    intersection.distance = std::numeric_limits<float>::max();
    intersection.triangleIndex = -1;
    intersection.sphereIndex = -1;

    bool intersects = false;
    float t0, t1; // distance
    glm::vec3 x0, x1;  // intersection points

    for(int t_i=0; t_i<triangles.size(); t_i++)
    {
        // Backface culling
        if(!isVisible(triangles[t_i].normal, dir))
            continue;

        if(triangles[t_i].intersect(start, dir, x0, t0) && 
           t0<intersection.distance)
        {   // valid intersection
            intersects = true;
            intersection.position = x0;
            intersection.distance = t0;
            intersection.triangleIndex = t_i;
        }
    }

    // Check intersection with spheres
    for(int s_i=0; s_i<spheres.size(); s_i++)
    {
        if(spheres[s_i].intersect(start, dir, x0, x1, t0, t1) &&
           ((t0>0.001 && t0<intersection.distance) || 
            (t1>0.001 && t1<intersection.distance) ))
        {
            intersects = true;
            intersection.position = t0>0? x0 : x1;
            intersection.distance = t0>0? t0 : t1;
            intersection.sphereIndex = s_i;
        }
    }
    return intersects;
}

void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result)
{
    int N = result.size();
    glm::vec2 step = glm::vec2(b-a) / float(max(N-1,1));
    glm::vec2 current( a );

    for(int i=0; i<N; i++) 
    {
        result[i] = current;
        current += step;
    }
}

float random_num(float min=-1.0, float max=1.0)
{
	float r = (float)rand() / (float)RAND_MAX;
	return min + r * (max - min);
}

glm::vec3 random_direction(const glm::vec3& n) 
{
    float z = sqrt(random_num(0));
    float r = sqrt(1.0 - z * z);
    float phi = 2.0 * PI * random_num(0);
    float x = r * cos(phi);
    float y = r * sin(phi);

    // local orthogonal coordinate system around n
    glm::vec3 w, u, v;
    w = n;
    u = glm::normalize(glm::cross(abs(w.x)>.1 ? glm::vec3(0,1,0) : 
                       glm::vec3(1,0,0), w));
    v = glm::cross(w, u);

    return x*u + y*v + z*w;
}

#endif
