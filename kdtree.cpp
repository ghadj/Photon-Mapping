#include "kdtree.h"

// Returns the square distance between points a and b
float KdTree::getSqrDist(const vec3 a, const vec3 b)
{
    float sqr_dist = 0.0;
    for(int i=0; i<3; i++)
        sqr_dist += (a[i] - b[i]) * (a[i] - b[i]);
    return sqr_dist;
}

void KdTree::buildNode(int* indices, int n_photons, int depth) 
{
    if (n_photons <= 0) return;

    // choose separation axis
    const int axis = depth % 3;

    // sort indices by coordination in the separation axis.
    sort(indices, indices + n_photons, [&](const int idx1, const int idx2) {
      return photons[idx1].position[axis] < photons[idx2].position[axis];
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

void KdTree::searchKNearestNode(int nodeIdx, const vec3& queryPoint, int k,
                                KNNQueue& queue) const 
{
    if (nodeIdx == -1 || nodeIdx >= nodes.size()) return;

    const Node& node = nodes[nodeIdx];

    // median point
    const Photon& median = photons[node.idx];

    // push to queue
    queue.emplace(getSqrDist(queryPoint, median.position), node.idx);

    // if size of queue is larger than k, pop queue
    if (queue.size() > k) 
    {
        queue.pop();
    }

    // if query point is lower than median, search left child
    // else, search right child
    const bool isLower = queryPoint[node.axis] < median.position[node.axis];
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
    const float dist_to_siblings = median.position[node.axis] - 
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

KdTree::KdTree(const Photon* photons, int nPhotons) 
{
    this->photons = photons;
    this->nPhotons = nPhotons;
}

void KdTree::buildTree() 
{
    // setup indices of photons
    vector<int> indices(nPhotons);
    iota(indices.begin(), indices.end(), 0);

    // build tree recursively
    buildNode(indices.data(), nPhotons, 0);
}

vector<int> searchKNearest(const vec3& queryPoint, int k, float& maxDist2)
    const 
{
    KNNQueue queue;
    searchKNearestNode(0, queryPoint, k, queue);

    vector<int> ret(queue.size());
    maxDist2 = 0;
    for (int i = 0; i < ret.size(); ++i) 
    {
        const auto& p = queue.top();
        ret[i] = p.second;
        maxDist2 = max(maxDist2, p.first);
        queue.pop();
    }
    return ret;
}
