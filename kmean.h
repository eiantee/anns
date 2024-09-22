//
// Created by caiqibang 00906710 on 2024/9/22.
//

#ifndef KMEAN_H
#define KMEAN_H
#include <limits>

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
};

template<int D>
void doKmean(cluster<D> & cluster, result<D> & result) {


}

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

#endif //KMEAN_H
