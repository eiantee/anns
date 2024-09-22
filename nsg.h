//
// Created by caiqibang 00906710 on 2024/9/21.
//

#ifndef NSG_H
#define NSG_H
#include <unordered_map>

template<int D, int M>
struct n_node {
    float value[D];
    int n_size;
    int neighbors[M];
    // center nearest
    n_node* parent;
};

template<int C>
struct nsg_candidate {
    int size;
    nsg_candidate() {
        size = 0;
    }
};


template<int D_SIZE, int D, int M, int C>
class Nsg {
public:
    n_node* root;
    n_node<D, M> nodes[D_SIZE];
    nsg_candidate<C> candidate[D_SIZE];

    void compute() {
        std::unordered_map<int, float> visited;
        for (int i=0; i<D_SIZE; i++) {
            visited.clear();
            for (int j = 0; j < nodes[i].parent->n_size; j++) {


        }
    }

    void searchOne(int start, float* query, std::unordered_map<int, float>& visited) {

    }
};

#endif //NSG_H
