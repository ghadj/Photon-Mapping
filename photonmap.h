#ifndef PHOTONMAP_H
#define PHOTONMAP_H 

#include "photon.h"
#include "kdtree.h"

using namespace std;

class PhotonMap
{
private:
    vector<Photon> photons;
    KdTree kdtree;

public:
    PhotonMap();
    void buildKdTree();
    void addPhoton();
};

#endif
