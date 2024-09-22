//
// Created by caiqibang 00906710 on 2024/9/21.
//

#ifndef NSG_H
#define NSG_H
#include <unordered_map>
#include <unordered_set>

#include "kmean.h"

template<int D, int M>
struct n_node {
    float value[D];
    int n_size;
    int neighbors[M];
    // center nearest
    n_node* parent;
    void addNeighbor(int i) {

    }
};

template<int C>
struct nsg_candidate {
    int size;
    int candi[C]{};
    float dis[C]{};
    nsg_candidate() {
        size = 0;
    }
};


template<int D_SIZE, int D, int M, int C>
class Nsg {
public:
    int size{D_SIZE};
    Nsg* upper;
    n_node* root;
    n_node<D, M> nodes[D_SIZE];
    nsg_candidate<C> candidate[D_SIZE];

    void compute() {
        std::unordered_map<int, float> dis_map;
        for (int i=0; i<D_SIZE; i++) {
            // TODO
            // prune
            assert(candidate[i].size > 0);
            nodes[i].addNeighbor(candidate[i].candi[0]);
            dis_map.emplace(candidate[i].candi[0], candidate[i].dis[0]);
            for (int j=1; j < candidate[i].size; j++) {
                dis_map.emplace(candidate[i].candi[j], candidate[i].dis[j]);
                for (int k=0; k<nodes[i].n_size; k++) {
                    // TODO
                    float temp_dis = vector_dis<D>(nodes[i].neighbors[k], candidate[i].candi[j]);
                    if (candidate[i].dis[j] > temp_dis && candidate[i].dis[j] > dis_map[nodes[i].neighbors[k]]) continue;
                    nodes[i].addNeighbor(candidate[i].candi[j]);
                }
            }
        }

        // tree
        std::unordered_set<int> to_travel;
        for (int i=0; i<upper->size; i++) {
            to_travel.clear();
            for (int)
            for (int k=0; k<nodes[i].size; k++) {
            

        }
    }

    void searchOne(int start, float* query, std::unordered_map<int, float>& visited) {

    }
};

#endif //NSG_H
