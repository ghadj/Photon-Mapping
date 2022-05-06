#ifndef KDTREE_H
#define KDTREE_H

#include <algorithm>
#include <glm/glm.hpp>
#include <queue>
#include "photon.h"

using namespace std;
using glm::vec3;

struct NeighborPhoton
{
    int index;  // index of photon
    float dist; // distance from query point
};

class KdTree 
{
private:
    struct Node 
    {
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

    // Returns the square distance between points a and b
    static float getSqrDist(const vec3 a, const vec3 b)
    {
        float sqr_dist = 0.0;
        for(int i=0; i<3; i++)
            sqr_dist += (a[i] - b[i]) * (a[i] - b[i]);
        return sqr_dist;
    }

    void buildNode(int* indices, int n_photons, int depth) 
    {
        if (n_photons <= 0) return;

        // choose separation axis
        const int axis = depth % 3;

        // sort indices by coordination in the separation axis.
        std::sort(indices, indices + n_photons, [&](const int idx1, const int idx2) {
          return photons[idx1].getDestination()[axis] < photons[idx2].getDestination()[axis];
        });

        // index of middle element of indices
        const int mid = (n_photons - 1) / 2;

        // add node to node array, remember index of current node(parent node)
        const int parentIdx = nodes.size();
        Node node;
        node.axis = axis;
        node.idx = indices[mid];
        nodes.push_back(node);

        // add left children to node array
        const int leftChildIdx = nodes.size();
        buildNode(indices, mid, depth + 1);

        // set index of left child on parent node
        // if size of nodes doesn't change, it means there is no left children
        if (leftChildIdx == nodes.size()) 
        {
            nodes[parentIdx].leftChildIdx = -1;
        } 
        else 
        {
            nodes[parentIdx].leftChildIdx = leftChildIdx;
        }

        // add right children to node array
        const int rightChildIdx = nodes.size();
        buildNode(indices + mid + 1, n_photons - mid - 1, depth + 1);

        // set index of right child on parent node
        // if size of nodes doesn't change, it means there is no right children
        if (rightChildIdx == nodes.size()) 
        {
            nodes[parentIdx].rightChildIdx = -1;
        } 
        else 
        {
            nodes[parentIdx].rightChildIdx = rightChildIdx;
        }
    }

    void searchKNearestNode(int nodeIdx, const vec3& queryPoint, int k,
                                    KNNQueue& queue) const 
    {
        if (nodeIdx == -1 || nodeIdx >= nodes.size()) return;

        const Node& node = nodes[nodeIdx];

        // median point
        const Photon& median = photons[node.idx];

        // push to queue
        queue.emplace(getSqrDist(queryPoint, median.getDestination()), node.idx);

        // if size of queue is larger than k, pop queue
        if (queue.size() > k) 
        {
            queue.pop();
        }

        // if query point is lower than median, search left child
        // else, search right child
        const bool isLower = queryPoint[node.axis] < median.getDestination()[node.axis];
        if (isLower) 
        {
            searchKNearestNode(node.leftChildIdx, queryPoint, k, queue);
        } 
        else 
        {
            searchKNearestNode(node.rightChildIdx, queryPoint, k, queue);
        }

        // at leaf node, if size of queue is smaller than k, or queue's largest
        // minimum distance overlaps sibblings region, then search siblings
        const float dist_to_siblings = median.getDestination()[node.axis] - 
                                       queryPoint[node.axis];

        if (queue.top().first > dist_to_siblings * dist_to_siblings) 
        {
            if (isLower) 
            {
                searchKNearestNode(node.rightChildIdx, queryPoint, k, queue);
            } 
            else 
            {
                searchKNearestNode(node.leftChildIdx, queryPoint, k, queue);
            }
        }
    }

public:
    KdTree(){}

    void setPhotons(const Photon* photons, int nPhotons) 
    {
        this->photons = photons;
        this->nPhotons = nPhotons;
    }

    void buildTree() 
    {
        // setup indices of photons
        vector<int> indices(nPhotons);
        iota(indices.begin(), indices.end(), 0);

        // build tree recursively
        buildNode(indices.data(), nPhotons, 0);
    }

    vector<NeighborPhoton> searchKNearest(const vec3& queryPoint, int k, 
                                          float& maxDist2) const 
    {
        KNNQueue queue;
        searchKNearestNode(0, queryPoint, k, queue);

        vector<NeighborPhoton> ret(queue.size());
        maxDist2 = 0;
        for (int i = 0; i < ret.size(); ++i) 
        {
            const auto& p = queue.top();
            ret[i].index = p.second;
            ret[i].dist = sqrt(p.first);
            maxDist2 = max(maxDist2, p.first);
            queue.pop();
        }
        return ret;
    }

};

#endif
