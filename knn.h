//
// Created by caiqibang 00906710 on 2024/9/21.
//

#ifndef KNN_H
#define KNN_H
#include <queue>
#include <unordered_set>
#include <bits/random.h>

#include "graph.h"

// M 近邻数量， D_SIZE 数据集大小
#define MAX_ITER 20
#define MAX_CANDI 1000

template<int M, int B>
struct knn_node {
    int n_size;
    int b_size;
    knn_node* next;
    int neighbor[M]{};
    int b_[M * B]{};
    knn_node() {
        n_size = 0;
        b_size = 0;
        next = nullptr;
    }
};

template<int C>
struct knn_candidate {
    int new_size;
    int old_size;
    knn_candidate* next;
    int candi_new[C]{};
    float dis_new[C]{};
    int candi_old[C]{};
    float dis_old[C]{};
    knn_candidate() {
        new_size = 0;
        old_size = 0;
        next = nullptr;
    }
};

template<int D_SIZE, int M, int C>
class Knn {
    public:
    knn_node<M, 8> node[D_SIZE];
    knn_candidate<C> candidate[D_SIZE];

    void compute() {
        // 随机
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, D_SIZE - 1);

        // 初始化随机neighbor
        std::unordered_set<int> randomNumbers;
        for (int i = 0; i < D_SIZE; i++) {
            while (randomNumbers.size() < M) {
                int num = dis(gen);
                if (num != i) {
                    randomNumbers.insert(num);
                }
            }
            node[i].add_neighbor(randomNumbers);
            node[i].addB(randomNumbers);
        }

        std::priority_queue<std::pair<float, int>> candidate_queue;
        std::unordered_set<int> b_;
        std::unordered_set<int> candidate_set;
        std::unordered_set<int> neighbor_set;
        int change_count = 0;
        for (int it = 0; it < MAX_ITER; it++) {
            change_count = 0;
            // 计算B_
            for (int i = 0; i < D_SIZE; i++) {
                node[i].setBSize(0);
            }
            for (int i = 0; i < D_SIZE; i++) {
                for (int j = 0; j < node[i].n_size; j++) {
                    node[node[i].getNeighbor(j)].addB(i);
                }
            }

            for (int i = 0; i < D_SIZE; i++) {
                b_.clear();
                candidate_set.clear();
                while (!candidate_queue.empty()) {
                    candidate_queue.pop();
                }
                for (int x = 0; x < node[i].getBSize(); x++) {
                    b_.insert(node[i].getB(x));
                }

                for (int b : b_) {
                    for (int x = 0; x < node[b].getBSize(); x++) {
                        candidate_set.insert(node[b].getB(x));
                    }
                }

                candidate_set.erase(i);
                for (int x = 0; x < candidate[i].getOldSize(); x++) {
                    // 去除上回已经对比过的
                    candidate_set.erase(candidate[i].getOld(x));
                    candidate_queue.emplace(candidate[i].getOldDis(x), candidate[i].getOld(x));
                }

                for (int c : candidate_set) {
                    float dis = vector_dis(i, dis);
                    if (candidate_queue.size() > MAX_CANDI && dis > candidate_queue.top().first) {
                        continue;
                    }
                    candidate_queue.emplace(dis, c);
                    if (candidate_queue.size() > MAX_CANDI) {
                        candidate_queue.pop();
                    }
                }

                // 更新neighbor 和candidate
                neighbor_set.clear();
                for (int x = 0; x < node[i].getNSize(); x++) {
                    neighbor_set.insert(node[i].getNeighbor(x));
                }
                node[i].setNSize(0);
                candidate[i].setOldSize(0);
                while (candidate_queue.size() > 0) {
                    candidate[i].addOld(candidate_queue.top());
                    if (candidate_queue.size() <= M) {
                        if (neighbor_set.find(candidate_queue.top().second) == neighbor_set.end()) {
                            change_count++;
                        }
                        node[i].addNeighbor(candidate_queue.top().second);
                    }
                    candidate_queue.pop();
                }
            }






        }



    }
};






#endif //KNN_H
