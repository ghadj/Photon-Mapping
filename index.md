## The Photon-Mapping Algorithm

The Photon-Maps is a global illumination algorithm and consists of two steps.

1. Building the photon map: emitting photons from the light sources into the
   scene and storing them in a photon map when they hit non-specular objects. 

2. The rendering pass: using the created photon map to extract information about
   incoming flux and reflected radiance at any point in the scene.

Purpose of this project is the implementation of the Photon-Mapping Algorithm on
a CPU renderer and test on a variation of the Cornell Box scene.

### Table of Contents

1. [The test environment](#1-the-test-environment)

2. [Photon tracing](#2-photon-tracing)
        
    2.1 [Light source & Photon Emission](#21-light-source-&-photon-emission)

    2.2 [Intersection with triangles and spheres](#22-intersection-with-triangles-and-spheres)

    2.3 [Reflection, refraction and absorption](#23-reflection-refraction-and-absorption)

3. [Photon-map](#3-photon-map)

    3.1 [The balanced kd-tree data structure](#31-the-balanced-kd-tree-data-structure)

4. [Rendering Pass](#4-rendering-pass)

    4.1 [Ray tracing](#41-ray-tracing) 
        
    4.2 [Estimation of radiance](#42-estimation-of-radiance)

5. [Discussion & Conclusion](#5-discussion-&-conclusion)
    
    5.1 [Evaluation](#51-evaluation)

    5.2 [Optimizations](#52-optimizations)


#### 1. The test environment

The implementation is tested on a variation of the Cornell Box. More
specifically the scene includes a diffuse point of light, a box and a glass
sphere. The glass sphere is added, in order to demonstrate the reflection
and refraction of light. For simplicity the rest surfaces are assumed to be
ideally diffused.

All the surfaces except of the sphere are represented as triangles, using three
coordinates, a normal vector and a color. The sphere is represented by the
coordinates of its center and the length of its radius.

#### 2. Photon Tracing

##### 2.1 Light source & Photon Emission 

During the photon tracing pass photons are emitted from the light sources of
the scene and according to their distribution of emissive power. For instance,
for a diffuse point of light, the distribution is uniform across random
directions of the source. For the purposes of this project, a
diffuse point of light is used, whose implementation is based on the pseudocode
of (Jensen & Christensen, 2000).

```
emit_photons_from_diffuse_point_light() {
    n_e = 0  //  number of emitted photons
    while (not enough photons) {
       // use simple rejection sampling to find 
       // diffuse photon direction
       do { 
            x = random number between -1 and 1
            y = random number between -1 and 1
            z = random number between -1 and 1
        } while ( x^2 + y^2 + z^2 > 1 )
        d = < x, y, z >
        p = light source position
        trace photon from p in direction d
        n_e = n_e + 1
    }
    scale power of stored photons with 1/n_e
}
```

##### 2.2 Intersection with triangles and spheres

Photons emitted from the light source intersect with the surfaces in the scene.

Assuming v_0, v_1, v_2 the coordinates of the edges of a triangle, in order to
describe a point in the triangle we construct a coordinate system that is
aligned with the triangle, with v0 to be the origin. Thus,

```
e_1 = v_1 - v_0
e_2 = v_2 - v_0
```

Any point r in the plane of the triangle can be described by u, v such that

```
r = v_0 + u*e_1 + v*e_2 (6)
```

For the points included within the triangle the following equalities apply:

```
0 < u
0 < v
u + v < 1
```

For each photon we store its origin point `o`  and direction `d`. Assuming `t` a
scalar coordinate describing the position on a photon's path, all the points r on the
path can be written as

```
r = o + td (10)
```

In order to find the intersection between the plane of the triangle and the line of the photon we combine the equations (6) and (10) and solve for the coordinates t, u, v:

```
v_0 + ue_1 + ve_2 = s + td (12)
```

The above equation can be formulated in 3x3 matrix notation as 

```
x = A^-1 b
```

where

```
x = (t u v )^T
A = (-d e_1 e_2)
b = s - v_0
```

Checking the inequalities (TODO) with the coordinates t, u, v for the intersection
point, we determine if the intersection is within the triangle.

In the case of intersection with sphere, a similar analytic solution was used. A
sphere can be described as 

```
x^2 + y^2 + z^2 = R^2
```

where x, y, z are the coordinates of a point P and R is the radius of the sphere. The above equation can be rewritten as 

```
P^2 - R^2 = 0
```

By substituting the equation of the photon's path with P, we get

```
|o + td|^2 - R^2 = 0
```

which developing the equation becomes

```
o^2 + (dt)^2 + 2odt - R^2 = o^2 + d^2t^2 + 2odt - R^2 = 0
```

The above equation is a quadratic function `f(x) = ax^2 + bx + c`, where with
`a=d^2`, `b=2od` and `c=o^2−R^2`. Solving for `f(x)=0` we can find the two
intersection points if exist.

```
x=\dfrac{-b\pm\sqrt{b^2-4ac}}{2a}\\
\Delta=b^2-4ac
```

Based on the discriminant (Δ), as illustrated in the figure below there
could be five different cases (from left to right):

1. Δ>0, there are two solutions, thus the path of the photon intersects at two
   points with the sphere. In addition, the distances of the intersection points
   from the path's origin are both positive, which means that the origin is
   outside the sphere and the direction of the photon is towards the sphere.

2. Δ=0, there is one solution, thus the path of the photon intersects at a
   single point with the sphere.

3. Δ>0, there are two solutions, thus the path of the photon intersects at two
   points with the sphere. However, in contrast with case (1), the distance of
   one of the intersection points from the photon path's origin is negative and
   the other is positive, which means that the origin is inside the sphere.

4. Δ<0, there are not any solutions, i.e., the path of the photon does not
   intersect with the sphere. 

5. Δ>0, there are two solutions, thus the path of the photon intersects at two
   points with the sphere. However, both intersection points have negative
   distance from the photon path's origin, which means that the origin is
   outside of the sphere and its direction points away from it.

For cases (4) and (5) the photon does not intersect with the
sphere.

[t0 and t1 represent the distance from the photon's path origin. Adopted from (Scratchapixel, 2014)](https://www.scratchapixel.com/images/upload/ray-simple-shapes/rayspherecases.png?)


##### 2.3 Reflection, refraction and absorption 

Once a photon is emitted from the light source, it is traced in the scene. Each
photon is stored at diffuse surfaces and can either be reflected, refracted by
the glass sphere or absorbed. In the case of diffuse surfaces, a simplification
of the statistical technique of *Russian roulette* (Avro & Kirk, 1990) is
used, by randomly deciding if a ray is absorbed or reflected. For instance, when
a photon hits a surface, other than the sphere, the photon is reflected if it
has not bounced in any surface so far, or if a random sample is less than the
predefined reflectance probability. Otherwise the photon is absorbed and no
longer traced.

```
if (photon.bounces == 0 || GetRandomNum(0) < 0.5) // Diffuse or absorb
    TracePhoton(p);
```

For the purposes of this project, we assume ideally diffused surfaces,
i.e., every ray of light incident on the surface is scattered equally in all
directions. The reflection of light is simulated using Monte Carlo integration,
by sampling the Bidirectional Reflectance Distribution Function (BRDF), which in
the case of ideally diffused surfaces, BRDF is uniform. The following pseudocode
is used to calculate a random new direction given the normal (n) of the surface
(Zhao, 2017).

```
uniform_random_dir(n) {
    z = sqrt(rand())
    r = sqrt(1.0 - z * z)
    φ = 2.0 * π * rand()
    x = r * cos(φ)
    y = r * sin(φ)
    [u, v, w] = create local (orthogonal) coordinate
                system around n
    return x*u + y*v + z*w
}
```

In the case a photon hits the sphere, then the Fresnel coefficient is used to
decide if the photon will be reflected or refracted. The Fresnel equation
describes how much light is reflected and how much is refracted, by considering
the angle of incidence. For example, when a photon hits the edges of the sphere,
it is most probable to be reflected, whereas incident photons in the center to
be refracted. Schlick's approximation formula was used to approximate the
contribution of the Fresnel factor (Schlick, 1994). The specular reflection
coefficient was calculated by:

[](add https://wikimedia.org/api/rest_v1/media/math/render/svg/7e790ea97c486d0d2d61dbca68be65c96aafd765)

Considering the approximated Fresnel coefficient R, a photon is randomly either
secularly reflected or refracted.

```
if(GetRandomNum() < R)
    Reflect photon
else
    Refract photon
```

The direction after specular reflection is calculated by 

```
Photon_Direction - 2 * (Surface_Normal . Photon_Direction) * Surface_Normal
```

Direction after refraction is calculated by 

$$ T = \dfrac{\eta_1}{\eta_2}(I + \cos(\theta_1)N) - N\sqrt{1 - \left( \dfrac{\eta_1}{\eta_2} \right) ^2 \sin^2(\theta_1)} $$

where theta_i is the angle of incidence, N the surface normal, n1 and n2 are the
indices of refraction which characterize the travel speed of light in different
materials. For this project we assumed index of refraction equal to 1.75, i.e.,
1.0 for the air and 1.75 for the glass materials.

The above equation can be rewritten as (Scratchapixel, 2014)

\begin{array}{l}
\eta = \dfrac{\eta_1}{\eta_2},\\
c_1 = \cos(\theta_1) = N \cdot I,\\
c_2 = \sqrt{1 - \left( \dfrac{n_1}{n_2} \right) ^2 \sin^2(\theta_1)} \rightarrow \sqrt{1 - \left( \dfrac{n_1}{n_2} \right) ^2 (1 - \cos^2(\theta_1))}
\end{array}

Thus, 

\begin{array}{l}
T = \eta(I + c_1 N) - N c_2,\\
T = \eta I + (\eta c_1 - c_2) N.
\end{array}

The calculation of refraction of light is calculated two times, when the photon
enters the sphere and when exits. Thus, in order to calculate the final
intersection point of the photon after intersecting with the sphere, first we
calculate the first direction of the photon after refraction, then given that
direction and the intersection point, a second intersection is calculated, where
the photon will exit, resulting in a new direction using the same above
equations. At the end, using the exit point and the final direction, we
calculate the closest intersection with other objects in the scene.


(Fresnel & calculation of power of new photon)

(picture of photons paths + picture of photons color at intersection points)

Caustics 

(picture of photons on caustics)


Considering all the above, photon tracing could be summarized in the following
pseudocode:

```
tracePhoton(p)
    find closest intersection

    if p intersects with the sphere
        if(GetRandomNum() < Fresnel Coefficient)
            Reflect photon and calculate new intersection point
        else
            Refract photon and calculate new intersection point

        add photon to the photon-map
    else
        if photon p bounces > 0
            add photon to the photon-map


    if photon p  bounces==0 OR GetRandomNum() < 0.5
        create photon p_new with origin from current intersection point and random direction
        tracePhoton(p_new)
```

#### 3. Photon-map

Photons are added to a list (a `vector` data structure) when they intersect with
diffuse surfaces either after they were refracted or reflected by the sphere or
other surface.

After building the photon map, we can later use it in the rendering pass to 
compute estimates of the incoming flux and the reflected radiance at
many points in the scene. To do this, it is necessary to locate the nearest
photons in the photon map. Because this operation is required to be done
very often, we can optimized the representation of the photon map before the
rendering pass, by using the *kd-tree* data structure.

##### 3.1 The balanced kd-tree data structure

A generic nearest neighbors search algorithm begins at the root of the kd-
tree, and adds photons to a list if they are within a certain distance.

The implementation of the kd-tree was adopted from (YumcyaWiz, 2021).

#### 4. Rendering Pass

Note that the photon map is view independent, and therefore it can be 
utilized to render the scene from any desired view. 

##### 4.1 Ray tracing 
        
##### 4.2 Estimation of radiance
      
(sphere vs triangles & direct lighting)

(pictures of different number of emitted photons)

(pictures of using photon-maps with and without direct lighting)

(pictures of different values of k-NN photons)

(picture with and without Fresnel coefficient)


### 5. Discussion & Conclusion

#### 5.1 Evaluation

#### 5.2 Optimizations

Optimizations: projection maps

### References

Avro, J., & Kirk, D. B. (1990). Particle transport and image synthesis in computer graphics. In SIGGRAPH’90: Proceedings of 17th Conference on Computer Graphics and Interactive Techniques (pp. 63-66).

Jensen, H. W., & Christensen, N. J. (2000). A practical guide to global illumination using photon maps. SIGGRAPH 2000 Course.

Schlick, C. (1994). An inexpensive BRDF model for physically‐based rendering. In Computer graphics forum (Vol. 13, No. 3, pp. 233-246). Edinburgh, UK: Blackwell Science Ltd.

Scratchapixel, (2014). Introduction to shading (reflection, refraction and Fresnel). Retrieved May 1, 2022, from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel 

Scratchapixel, (2014). A minimal ray-tracer: Rendering simple shapes (ray-sphere intersection). Retrieved May 1, 2022, from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection 

YumcyaWiz, (2021) Photon_mapping Retrieved May 1, 2022 from Github repository https://github.com/yumcyaWiz/photon_mapping

Zhao, S. (2017). Rendering Equation & Monte Carlo Path Tracing. CS295, Spring 2017, University of California, Irvine. 

### Other useful links

+ [CMU Lecture slides on Photon mapping](https://www.cs.cmu.edu/afs/cs/academic/class/15462-s12/www/lec_slides/lec18.pdf)

