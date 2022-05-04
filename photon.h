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
        : normal(normal), position(position), energy(energy), bounces(0){}

    vec3 getNormal() const
    {
        return this->normal;
    }

    vec3 getPosition() const
    {
        return this->position;
    }

    int getBounces() const
    {
        return this->bounces;
    }

    //TODO
    //float operator[](int i) const { return position[i]; }
};

#endif
