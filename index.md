## The Photon-Mapping Algorithm

Purpose of this project is the implementation of the Photon-Mapping Algorithm on
a CPU renderer and test on a variation of the Cornell Box scene.

### Table of Contents

1. [Implementation](#implementation)

    1.1 [Photon Tracing](#photon-tracing)

        1.1.1 [The balanced kd-tree data structure](#the-balanced-kd-tree-data-structure)

    1.2 [Rendering Pass](#rendering-pass)

2. [Evaluation](#evaluation)

3. [Discussion & Conclusion](#discussion-&-conclusion)

### 1. Implementation

The Photon-Maps is a global illumination algorithm and consists of two steps.

1. Building the photon map: emitting photons from the light sources into the
   scene and storing them in a photon map when they hit non-specular objects. 

2. The rendering pass: using the created photon map to extract information about
   incoming flux and reflected radiance at any point in the scene.

#### 1.1 Photon Tracing

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
        do { // use simple rejection sampling to find diffuse photon direction
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

Each photon is stored (at diffuse surfaces) and absorbed(A) or reflected(R) or
transmitted(T). The following equality applies:

    A+R+T = 1.0

A simple approach to this, would be given an incoming photon with power Φ, to
reflect the photon with power R\*Φ and transmit is with power T\*Φ.  However,
there is the risk to end up having too many low-powered photons, which would be
a waste of resources. The statistical technique of *Russian roulette* (Avro &
Kirk, 1990) can be used to solve this issue, by randomly deciding if a ray is
absorbed or reflected. For example, given the surface reflectance R,
transmittance T and incoming photon power Φ:

```
r = random();
if (r < R)
    reflect photon with power Φ
else if (r < T)
    transmit photon with power Φ
else
    photon is absorbed
```

That way we keep fewer photons with similar power and terminate un-important
ones.

Note that the photon map is view independent, and therefore it can be 
utilized to render the scene from any desired view. 

After building the photon map, we can later use it in the rendering pass to 
compute estimates of the incoming flux and the reflected radiance at
many points in the scene. To do this, it is necessary to locate the nearest
photons in the photon map. Because this operation is required to be done
very often, we can optimized the representation of the photon map before the
rendering pass, by using the *kd-tree* data structure.

Thus, extending the pseudocode for emitting photons, we have:

```
while (we want more photons) {
    emit a photon
    while (photon hits a surface) {
        store photon
        use Russian Roulette to scatter photon
    }
}
build the kd-tree 
```

##### The balanced kd-tree data structure

A generic nearest neighbors search algorithm begins at the root of the kd-
tree, and adds photons to a list if they are within a certain distance.

#### 1.2 Rendering Pass

Ray tracer... 

### 2. Evaluation

### 3. Discussion & Conclusion

optimizations: projection maps

### References

Jensen, H. W., & Christensen, N. J. (2000). A practical guide to global illumination using photon maps. SIGGRAPH 2000 Course Notes CD-ROM.

Avro, J., & Kirk, D. B. (1990). Particle transport and image synthesis in computer graphics. In SIGGRAPH’90: Proc. of 17th Conf. on Computer Graphics and Interactive Techniques (pp. 63-66).

### Other useful links

https://www.cs.cmu.edu/afs/cs/academic/class/15462-s12/www/lec_slides/lec18.pdf
