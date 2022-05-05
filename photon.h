#ifndef PHOTON_H
#define PHOTON_H

#include <glm/glm.hpp>

using glm::vec3;

class Photon
{
private:
    vec3 normal;
    vec3 source;
    vec3 destination;
    vec3 energy;
    int bounces;
     
public:
    Photon(vec3 normal, vec3 source, vec3 energy, int bounces=0) 
        : normal(normal), source(source), energy(energy), bounces(bounces){}

    vec3 getNormal() const
    {
        return this->normal;
    }

    vec3 getSource() const
    {
        return this->source;
    }

    void setDestination(vec3 d)
    {
        this->destination = d;
    }

    vec3 getDestination() const
    {
        return this->destination;
    }

    vec3 getEnergy() const
    {
        return this->energy;
    }
    
    int getBounces() const
    {
        return this->bounces;
    }

    //TODO
    //float operator[](int i) const { return position[i]; }
};

#endif
