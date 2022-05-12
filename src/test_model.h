#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

#include "shapes.h"

using glm::vec3;

// Defines colors
const vec3 BLUE(0.35f, 0.35f, 0.95f);
const vec3 RED(0.95f, 0.35f, 0.35f);
const vec3 WHITE(0.95f, 0.95f, 0.95f);

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel(std::vector<Triangle> &triangles, std::vector<Sphere> &spheres)
{
	triangles.clear();
	triangles.reserve(5 * 2 * 3);

	// Define Room
	float L = 555; // Length of Cornell Box side.

	vec3 A(L, 0, 0);
	vec3 B(0, 0, 0);
	vec3 C(L, 0, L);
	vec3 D(0, 0, L);

	vec3 E(L, L, 0);
	vec3 F(0, L, 0);
	vec3 G(L, L, L);
	vec3 H(0, L, L);

	// Floor
	triangles.push_back(Triangle(C, B, A, WHITE));
	triangles.push_back(Triangle(C, D, B, WHITE));

	// Left wall
	triangles.push_back(Triangle(A, E, C, RED));
	triangles.push_back(Triangle(C, E, G, RED));

	// Right wall
	triangles.push_back(Triangle(F, B, D, BLUE));
	triangles.push_back(Triangle(H, F, D, BLUE));

	// Ceiling
	triangles.push_back(Triangle(E, F, G, WHITE));
	triangles.push_back(Triangle(F, H, G, WHITE));

	// Back wall
	triangles.push_back(Triangle(G, D, C, WHITE));
	triangles.push_back(Triangle(G, H, D, WHITE));

	// Front wall
	triangles.push_back(Triangle(F, E, B, WHITE));
	triangles.push_back(Triangle(B, E, A, WHITE));

	// Define the sphere
	spheres.push_back(Sphere(vec3(120, 100, 90), 100.0));

	/*/ Define the Short block
	A = vec3(290,0,114);
	B = vec3(130,0, 65);
	C = vec3(240,0,272);
	D = vec3( 82,0,225);

	E = vec3(290,165,114);
	F = vec3(130,165, 65);
	G = vec3(240,165,272);
	H = vec3( 82,165,225);

	// Front
	triangles.push_back(Triangle(E, B, A, WHITE));
	triangles.push_back(Triangle(E, F, B, WHITE));

	// Front
	triangles.push_back(Triangle(F, D, B, WHITE));
	triangles.push_back(Triangle(F, H, D, WHITE));

	// BACK
	triangles.push_back(Triangle(H, C, D, WHITE));
	triangles.push_back(Triangle(H, G, C, WHITE));

	// LEFT
	triangles.push_back(Triangle(G, E, C, WHITE));
	triangles.push_back(Triangle(E, A, C, WHITE));

	// TOP
	triangles.push_back(Triangle(G, F, E, WHITE));
	triangles.push_back(Triangle(G, H, F, WHITE));
	*/

	// Define the Tall block
	A = vec3(423, 0, 247);
	B = vec3(265, 0, 296);
	C = vec3(472, 0, 406);
	D = vec3(314, 0, 456);

	E = vec3(423, 330, 247);
	F = vec3(265, 330, 296);
	G = vec3(472, 330, 406);
	H = vec3(314, 330, 456);

	// Front
	triangles.push_back(Triangle(E, B, A, WHITE));
	triangles.push_back(Triangle(E, F, B, WHITE));

	// Front
	triangles.push_back(Triangle(F, D, B, WHITE));
	triangles.push_back(Triangle(F, H, D, WHITE));

	// BACK
	triangles.push_back(Triangle(H, C, D, WHITE));
	triangles.push_back(Triangle(H, G, C, WHITE));

	// LEFT
	triangles.push_back(Triangle(G, E, C, WHITE));
	triangles.push_back(Triangle(E, A, C, WHITE));

	// TOP
	triangles.push_back(Triangle(G, F, E, WHITE));
	triangles.push_back(Triangle(G, H, F, WHITE));

	// Scale triangles to the volume [-1,1]^3
	for (size_t i = 0; i < triangles.size(); ++i)
	{
		triangles[i].v0 *= 2 / L;
		triangles[i].v1 *= 2 / L;
		triangles[i].v2 *= 2 / L;

		triangles[i].v0 -= vec3(1, 1, 1);
		triangles[i].v1 -= vec3(1, 1, 1);
		triangles[i].v2 -= vec3(1, 1, 1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].ComputeNormal();
	}

	// Scale spheres to the volume [-1,1]^3
	for (size_t i = 0; i < spheres.size(); ++i)
	{
		spheres[i].radius *= 2 / L;
		spheres[i].center *= 2 / L;

		spheres[i].center -= vec3(1, 1, 1);
		spheres[i].center.x *= -1;
		spheres[i].center.y *= -1;
	}
}

#endif
