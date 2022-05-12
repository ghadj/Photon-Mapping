/**
 * @file shapes.h
 * @author Georgios Hadjiantonis
 * @brief Definitions of surface/object classes used in the scene.
 *
 */
#ifndef SHAPES_H
#define SHAPES_H

#include <glm/glm.hpp>

using glm::mat3;
using glm::vec3;

/**
 * @brief Class of a triangular surface.
 *
 */
class Triangle
{
public:
    vec3 v0;
    vec3 v1;
    vec3 v2;
    vec3 normal;
    vec3 color;

    Triangle(vec3 v0, vec3 v1, vec3 v2, vec3 color)
        : v0(v0), v1(v1), v2(v2), color(color)
    {
        ComputeNormal();
    }

    /**
     * @brief Computes the normal vector of the triangle.
     *
     */
    void ComputeNormal()
    {
        vec3 e1 = v1 - v0;
        vec3 e2 = v2 - v0;
        normal = glm::normalize(glm::cross(e2, e1));
    }

    /**
     * @brief Given an origin point and a direction, checks if the ray
     *        intersects with the triangle.
     *
     * @param start  Origin of the ray.
     * @param dir    Direction vector of the ray.
     * @param x0     Intersection point.
     * @param t0     Distance from origin to intersection point.
     * @return true  If the ray intersects with the triangle.
     * @return false If the ray does not intersect with the triangle.
     */
    bool Intersect(const vec3 start, const vec3 dir, vec3 &x0, float &t0) const
    {
        vec3 e1, e2, b, x;
        e1 = this->v1 - this->v0;
        e2 = this->v2 - this->v0;
        b = start - this->v0;
        mat3 A(-dir, e1, e2);
        x = glm::inverse(A) * b; // x:t, y:u, z:v

        // Add equallity so that it includes points on the edges of the triangle
        if (x.x >= 0 && x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1)
        {
            x0 = this->v0 + x.y * e1 + x.z * e2;
            t0 = x.x;
            return true;
        }
        return false;
    }
};

/**
 * @brief Class of a sphere object.
 *
 */
class Sphere
{
public:
    vec3 center;
    float radius;

    Sphere(vec3 center, float radius)
        : center(center), radius(radius) {}

    /**
     * @brief Given an origin point and a direction, checks if the ray
     *        intersects with the sphere.
     *
     * Note that in the case of two intersection points, and the one is behind
     * the origin of the ray, in regards to its direction, the intersection
     * point will correspond to parameters x0 and t0.
     *
     * @param start  Origin of the ray.
     * @param dir    Direction vector of the ray.
     * @param x0     Intersection point.
     * @param x1     Intersection point
     * @param t0     Distance from origin to intersection point x0.
     * @param t1     Distance from origin to intersection point x1.
     * @return true  If the ray intersects with the triangle.
     * @return false If the ray does not intersect with the triangle.
     */
    bool Intersect(const vec3 &start, const vec3 &dir,
                   vec3 &x0, vec3 &x1,         // intersection points
                   float &t0, float &t1) const // distance from origin
    {
        vec3 center_trn = start - this->center; // translated center of sphere
        float A = glm::dot(dir, dir);
        float B = 2.0 * glm::dot(dir, center_trn);
        float C = glm::dot(center_trn, center_trn) - pow(this->radius, 2.0);
        float Delta = pow(B, 2.0) - 4 * A * C;

        if (Delta == 0.0) // there is only one solution
        {
            t0 = t1 = -0.5 * B / A;
            x0 = start + t0 * dir;
            x1 = start + t1 * dir;
            return true;
        }
        else if (Delta > 0.0) // there are two solutions
        {
            t0 = (-B + sqrt(Delta)) / (2 * A);
            t1 = (-B - sqrt(Delta)) / (2 * A);

            if (t0 < 0 && t1 < 0)
                return false;

            // Note: now there is the possibility one of t0 or t1 to be negative.
            // In that case the origin of the ray is inside the sphere.

            // Swap so that t1 corresponds to the farthest intersection point
            if (t0 > t1)
                std::swap(t0, t1);

            x0 = start + t0 * dir;
            x1 = start + t1 * dir;

            return true;
        }
        return false; // no intersection
    }
};

#endif
