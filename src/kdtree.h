/**
 * @file kdtree.h
 * @author Georgios Hadjiantonis
 * @brief Definition of a kd-tree.
 * 
 */
#ifndef KDTREE_H
#define KDTREE_H

#include <algorithm>
#include <queue>
#include <numeric>

#include <glm/glm.hpp>

#include "photon_tracer.h"

using namespace std;
using glm::vec3;

/**
 * @brief Instance of neighbor photon returned by the search.
 * 
 */
struct NeighborPhoton
{
    int index;  // index of the corresponding photon in the given array.
    float dist; // distance from query point
};

/**
 * @brief Definition of a kd-tree.
 * 
 */
class KdTree
{
private:
    struct Node
    {
        int axis;            // separation axis(x=0, y=1, z=2)
        int idx;             // index of median point
        int left_child_idx;  // index of left child
        int right_child_idx; // index of right child

        Node() : axis(-1), idx(-1), left_child_idx(-1), right_child_idx(-1) {}
    };

    vector<Node> nodes;    // array of tree nodes
    int num_photons;       // number of photons
    const Photon *photons; // pointer to array of photons

    using KNNQueue = priority_queue<pair<float, int>>;

    /**
     * @brief Returns the square distance between points a and b.
     * 
     * @param a      Coordinates of point a.
     * @param b      Coordinates of point b.
     * @return float Square distance between a and b.
     */
    static float GetSqrDist(const vec3 a, const vec3 b)
    {
        float sqr_dist = 0.0;
        for (int i = 0; i < 3; i++)
            sqr_dist += (a[i] - b[i]) * (a[i] - b[i]);
        return sqr_dist;
    }

    /**
     * @brief Recursively builds each node of the tree.
     * 
     * @param indices   Indices of photons in the array.
     * @param n_photons Number of photons left.
     * @param depth     Current depth in the tree.
     */
    void BuildNode(int *indices, int n_photons, int depth)
    {
        if (n_photons <= 0)
            return;

        // Choose separation axis
        const int axis = depth % 3;

        // Sort indices by coordination in the separation axis.
        std::sort(indices, indices + n_photons, [&](const int idx1, const int idx2)
                  { return photons[idx1].destination[axis] < photons[idx2].destination[axis]; });

        // Index of middle element of indices
        const int mid = (n_photons - 1) / 2;

        // Add node to node array, remember index of current node(parent node)
        const int parent_idx = nodes.size();
        Node node;
        node.axis = axis;
        node.idx = indices[mid];
        nodes.push_back(node);

        // Add left children to node array
        const int left_child_idx = nodes.size();
        BuildNode(indices, mid, depth + 1);

        // Set index of left child on parent node
        // If size of nodes doesn't change, it means there is no left children
        if (left_child_idx == nodes.size())
        {
            nodes[parent_idx].left_child_idx = -1;
        }
        else
        {
            nodes[parent_idx].left_child_idx = left_child_idx;
        }

        // Add right children to node array
        const int right_child_idx = nodes.size();
        BuildNode(indices + mid + 1, n_photons - mid - 1, depth + 1);

        // Set index of right child on parent node
        // If size of nodes doesn't change, it means there is no right children
        if (right_child_idx == nodes.size())
        {
            nodes[parent_idx].right_child_idx = -1;
        }
        else
        {
            nodes[parent_idx].right_child_idx = right_child_idx;
        }
    }

    /**
     * @brief Given the query point returns the k nearest photons in a priority queue.
     *
     * @param node_idx    Index of current node. 
     * @param query_point Point of interest.
     * @param k           Number of photons to be returned.
     * @param queue       Priority queue based on the distance of photons from
     *                    questy point.
     */
    void SearchKNearestNode(int node_idx, const vec3 &query_point, int k,
                            KNNQueue &queue) const
    {
        if (node_idx == -1 || node_idx >= nodes.size())
            return;

        const Node &node = nodes[node_idx];

        // Median point
        const Photon &median = photons[node.idx];

        // Push to queue
        queue.emplace(GetSqrDist(query_point, median.destination), node.idx);

        // If size of queue is larger than k, pop queue
        if (queue.size() > k)
        {
            queue.pop();
        }

        // If query point is lower than median, search left child
        // else, search right child
        const bool is_lower = query_point[node.axis] < median.destination[node.axis];
        if (is_lower)
        {
            SearchKNearestNode(node.left_child_idx, query_point, k, queue);
        }
        else
        {
            SearchKNearestNode(node.right_child_idx, query_point, k, queue);
        }

        // At leaf node, if size of queue is smaller than k, or queue's largest
        // minimum distance overlaps siblings region, then search siblings
        const float dist_to_siblings = median.destination[node.axis] -
                                       query_point[node.axis];

        if (queue.top().first > dist_to_siblings * dist_to_siblings)
        {
            if (is_lower)
            {
                SearchKNearestNode(node.right_child_idx, query_point, k, queue);
            }
            else
            {
                SearchKNearestNode(node.left_child_idx, query_point, k, queue);
            }
        }
    }

public:
    KdTree() {}

    /**
     * @brief Sets pointer to array of photons.
     * 
     * @param photons     Pointer to array of photons.
     * @param num_photons Total number of photons.
     */
    void SetPhotons(const Photon *photons, int num_photons)
    {
        this->photons = photons;
        this->num_photons = num_photons;
    }

    /**
     * @brief Given the array of photons, recursively builds the kd-tree.
     * 
     */
    void BuildTree()
    {
        // Setup indices of photons
        vector<int> indices(num_photons);
        iota(indices.begin(), indices.end(), 0);

        // Build tree recursively
        BuildNode(indices.data(), num_photons, 0);
    }

    /**
     * @brief Given the query point returns the k nearest photons and the
     *        maximum squared distance from the farthest photon in the list.
     * 
     * @param query_point             Query point.
     * @param k                       Number of photons to be returned.
     * @param max_dist2               Max squared distance of the k-th photon.
     * @return vector<NeighborPhoton> Vector of the k nearest photons.
     */
    vector<NeighborPhoton> SearchKNearest(const vec3 &query_point, int k,
                                          float &max_dist2) const
    {
        KNNQueue queue;
        SearchKNearestNode(0, query_point, k, queue);

        vector<NeighborPhoton> ret(queue.size());
        max_dist2 = 0;
        for (int i = 0; i < ret.size(); ++i)
        {
            const auto &p = queue.top();
            ret[i].index = p.second;
            ret[i].dist = sqrt(p.first);
            max_dist2 = max(max_dist2, p.first);
            queue.pop();
        }
        return ret;
    }
};

#endif
