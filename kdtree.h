#ifndef KDTREE_H
#define KDTREE_H

#include <glm/glm.hpp>
#include <queue>
#include "photon.h"

using namespace std;
using glm::vec3;

class KdTree {
private:
    struct Node {
        int axis;           // separation axis(x=0, y=1, z=2)
        int idx;            // index of median point
        int leftChildIdx;   // index of left child
        int rightChildIdx;  // index of right child

        Node() : axis(-1), idx(-1), leftChildIdx(-1), rightChildIdx(-1) {}
    };

    vector<Node> nodes;     // array of tree nodes
    int nPhotons;           // number of photons
    const Photon* photons;  // pointer to array of photons

    using KNNQueue = priority_queue<pair<float, int>>;

    float getSqrDist(const vec3 a, const vec3 b);
    void buildNode(int* indices, int n_photons, int depth);
    void searchKNearestNode(int nodeIdx, const vec3& queryPoint, 
                            int k, KNNQueue& queue) const;
public:
    KdTree(const Photon* photons, int nPhotons);
    void buildTree(); 
    vector<int> searchKNearest(const vec3& queryPoint, int k, 
                               float& maxDist2) const;
};

#endif
