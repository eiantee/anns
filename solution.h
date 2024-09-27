//
// Created by caiqibang 00906710 on 2024/9/27.
//

#ifndef CPP_SOLUTION_H
#define CPP_SOLUTION_H
#include <vector>
#include "layer.h"
#include "kmean.h";

template<int D>
struct handler{
    float* base;
    VBase<D> * v_base;
    ClusterIndex * cluster_index;
    LayIndex * layer_index;
    std::vector<std::vector<cluster<D> *>>* thread_k_mean;
    Layer<D, 10>* layer0;
    Layer<D, 10>* layer1;
    Layer<D, 10>* layer2;
    Layer<D, 10>* layer_base;
};

namespace global {
    inline static int DIM;
    inline static handler<512>* handler_512;
    inline static handler<1024>* handler_1024;
    inline static handler<2048>* handler_2048;
    inline static handler<TEST_DIM>* handler_test;

    template<int D>
    inline static handler<D>* get_handler() {
        if (D == 512) {
            return handler_512;
        }
        if (D == 1024) {
            return handler_1024;
        }
        if (D == 2048) {
            return handler_2048;
        }
        if (D == TEST_DIM) {
            return handler_test;
        }
        return nullptr;
    }
}

class Solution {
public:
    void build(int d, const std::vector<float>& base);
    void search(const std::vector<float>& query, int* res);
};

#endif //SOLUTION_H
