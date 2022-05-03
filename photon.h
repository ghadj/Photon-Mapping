#ifndef PHOTON_H
#define PHOTON_H

#include <glm/glm.hpp>

using glm::vec3;

class Photon
{
private:
    vec3 normal;
    vec3 position;
    vec3 energy;
    int bounces;
     
public:
    Photon(vec3 normal, vec3 position, vec3 energy) 
        : normal(normal), position(position), energy(energy), bounces(1){}

    //TODO
    //float operator[](int i) const { return position[i]; }
};


#endif
