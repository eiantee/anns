//
// Created by caiqibang 00906710 on 2024/9/22.
//
#pragma once

#ifndef KMEAN_H
#define KMEAN_H
#include <limits>
#include <vector>

#include "config.h"

#define MAX_KMEAN_ITER 20
#define ITER_THRESHOLD 0.1

template<int D>
struct cluster {
	float * center_value{nullptr};
    int size;
    int nearest_id;
    float * base;
    int * id_index;
    bool root;
    int offset{0};
    int base_line;

    float* at(int id) {
        return base + base_line * getRealId(id) + offset;
    }

    inline int getRealId(int id) {
        if (root) {
            return id;
        }
        return id_index[id];
    }

    cluster() {
        size = 0;
        base = nullptr;
        id_index = nullptr;
    }

    cluster(bool root_in, int offset_in, int base_line_in, int * id_index_in, float* base_in, int size_in) {
        root = root_in;
        offset = offset_in;
        base_line = base_line_in;
        size = size_in;
        id_index = id_index_in;
        base = base_in;
    }

    ~cluster() {
        if (id_index != nullptr) {
            delete[] id_index;
        }
        if (center_value != nullptr) {
            delete[] center_value;
        }
    }
};

template<int D>
struct center {
    int id;
    int nearest;
    int count;
    float value[D];
};

template<int D>
struct result {
    int center_size;
    int cluster_size;
    center<D> *centers;
    int * index;
    result(int center_s, int cluster_s) {
        center_size = center_s;
        cluster_size = cluster_s;
        centers = new cluster<D>[center_size];
        index = new int[cluster_size];
    }

    void choose(int data_id, int center_id) {
        index[data_id] = center_id;
    }

    ~result() {
        delete[] centers;
        delete[] index;
    }

    std::vector<cluster<D>*> split(cluster<D>& parent) {
        std::vector<cluster<D>*> res;
        for (int i = 0; i < center_size; i++) {
            if (centers[i].count < 1) continue;
            int* id_index = new int[centers[i].count];
            int cc = 0;
            int near = 0;
            for (int k = 0; k < cluster_size; k++) {
                if (index[k] == i) {
                    if (k == centers[i].nearest) {
                        near = cc;
                    }
                    id_index[cc] = parent.getRealId(k);
                    cc++;
                }
            }
            cluster<D>* tt = new cluster<D>(false, parent.offset, parent.base_line, id_index,
                parent.base, centers[i].count);
            tt->nearest_id = near;
            float* cv = new float[D];
            for (int v = 0; v < D; v++) {
                cv[v] = centers[i].value[v];
            }
            tt->center_value = cv;
            res.push_back(tt);
        }
        return res;
    }
};

template<int D>
void choose(std::function<float(int, int)> distance, result<D> & result) {
    // 随机
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, result.cluster_size - 1);

    int ids[result.center_size];
    int times = log2(result.cluster_size);
    ids[0] = dis(gen);

    for (int i = 1; i < result.center_size; i++) {
        int max_id = 0;
        float max_dis = 0.0f;
        for (int t = 0; t < times; t++) {
            float min_dis = std::numeric_limits<float>::max();
            int temp_id = dis(gen);
            for (int j = 0; j < i; j++) {
                float temp_dis = distance(temp_id, j);
                if (temp_dis < min_dis) {
                    min_dis = temp_dis;
                }
            }
            if (min_dis > max_dis) {
                max_dis = min_dis;
                max_id = temp_id;
            }
        }
        ids[i] = max_id;
    }

    for (int i = 0; i < result.center_size; i++) {
        result.centers[i].id = ids[i];
    }
}

template<int D>
void doKmean(cluster<D>& cluster, result<D>& result) {
    choose<D>([&cluster](int a, int b) {return vector_dis<D>(cluster.at(a), cluster.at(b));}, result);
    // init 初始向量
    for (int i = 0; i < result.center_size; i++) {
        for (int j = 0; j < D; j++) {
            result.centers[i].value[j] = *(cluster.at(result.centers[i].id) + j);
        }
    }

    double temp_center[result.center_size][D];
    int temp_count[result.center_size];
    // 分配系数, 遇到count 0 要从多的分配出来
    int weight[result.center_size];
    int nearest_id[result.center_size];
    float nearest_dis[result.center_size];
    double diff = 0.0;

    for (int iter = 0; iter < MAX_KMEAN_ITER; iter++) {
        // init
        diff = 0.0;
        for (int i = 0; i < result.center_size; i++) {
            temp_count[i] = 0;
            weight[i] = 1;
            nearest_dis[i] = std::numeric_limits<float>::max();
            for (int j = 0; j < D; j++) {
                temp_center[i][j] = 0.0;
            }
        }

        // 选择和计算中心
        for (int i = 0; i < result.cluster_size; i++) {
            int center_index = 0;
            float min_dis = std::numeric_limits<float>::max();
            for (int j = 0; j < result.center_size; j++) {
                float temp_dis = vector_dis<D>(result.centers[j].value, cluster.at(i));
                if (temp_dis < min_dis) {
                    min_dis = temp_dis;
                    center_index = j;
                }
            }
            result.choose(i, center_index);
            if (min_dis < nearest_dis[center_index]) {
                nearest_dis[center_index] = min_dis;
                nearest_id[center_index] = i;
            }
            temp_count[center_index] = temp_count[center_index] + 1;
            addf<D>(temp_center[center_index], cluster.at(i));
        }

        // 重新计算center
        for (int j = 0; j < result.center_size; j++) {
            float temp[D];
            if (temp_count[j] < 1) continue;
            scaled<D>(temp_center[j], 1.0/temp_count[j]);
            for (int k = 0; k < D; k++) {
                temp[k] = temp_center[j][k];
            }

            diff += vector_dis<D>(temp, result.centers[j].value);
            for (int k = 0; k < D; k++) {
                result.centers[j].value[k] = temp[k];
            }
        }

        if (diff < ITER_THRESHOLD) {
            break;
        }
    }
    for (int j = 0; j < result.center_size; j++) {
        result.centers[j].nearest = nearest_id[j];
        result.centers[j].count = temp_count[j];
    }
}

#endif //KMEAN_H
