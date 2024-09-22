//
// Created by caiqibang 00906710 on 2024/9/22.
//

#ifndef KMEAN_H
#define KMEAN_H
#include <limits>

#define MAX_KMEAN_ITER 20
#define ITER_THRESHOLD 0.1

template<int D>
struct cluster {
    int size;
    cluster * parent;
    int * index;

    float* at(int id) {
        index[id] = id;
    }


};

template<int D>
struct center {
    int id;
    int nearest;
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
float vector_dis(float* a, float* b) {
}

template<int D>
void vector_add(double* a, float* b) {
}

template<int D>
void vector_scale(double* a, double b) {
}


template<int D>
void doKmean(cluster<D> & cluster, result<D> & result) {
    choose([](int a, int b) { return vector_dis<D>(cluster.at(a), cluster.at(b));}, result);
    // init 初始向量
    for (int i = 0; i < result.center_size; i++) {
        for (int j = 0; j < D; j++) {
            result.centers[i].value[j] = cluster.at(result.centers[i].id) + j;
        }
    }

    double temp_center[result.center_size][D];
    int temp_count[result.center_size];
    int nearest_id[result.center_size];
    float nearest_dis[result.center_size];

    for (int iter = 0; iter < MAX_KMEAN_ITER; iter++) {

        // init
        for (int i = 0; i < result.center_size; i++) {
            temp_count[i] = 0;
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
                float temp_dis = vector_dis(result.centers[j].value, cluster.at(i));
                if (temp_dis < min_dis) {
                    min_dis = temp_dis;
                    center_index = j;
                }
            }
            result.choose(i, center_index);
            if (min_dis < nearest_dis[i]) {
                nearest_dis[center_index] = min_dis;
                nearest_id[center_index] = i;
            }
            temp_count[center_index] = temp_count[center_index] + 1;
            vector_add<D>(temp_center[center_index], cluster.at(i));
        }

        // 重新计算center
        double diff = 0.0;
        for (int j = 0; j < result.center_size; j++) {
            vector_scale<D>(temp_center[j], 1.0 / temp_count[j]);
            float temp[D];
            for (int k = 0; k < D; k++) {
                temp[k] = temp_center[j][k];
            }
            diff += vector_dis<D>(temp, result.centers[j].value);
            for (int k = 0; k < D; k++) {
                result.centers[j].value[k] = temp[k];
            }
        }

        if (diff < ITER_THRESHOLD) break;
    }

    for (int j = 0; j < result.center_size; j++) {
        result.centers[j].nearest = nearest_id[j];
    }
}

#endif //KMEAN_H
