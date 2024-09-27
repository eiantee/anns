//
// Created by caiqibang 00906710 on 2024/9/26.
//

#ifndef LAYER_H
#define LAYER_H
#include <array>
#include <atomic>
#include <queue>
#include <stdlib.h>
#include <unordered_set>
#include <vector>

#include "config.h"

#define  MAX_M 200
// K * OVER_LAP
#define  OVER_LAP 50
#define LAYER_DEP 5
#define V_MAX 100000

namespace layer {
    thread_local static short layer_thread_id = 0;
    static short init = -1;
}

template<int D>
struct VBase {
    float * v_b;
    std::atomic<int> size{0};
    VBase() {
        v_b = static_cast<float *>(aligned_alloc(CACHE_LINE, sizeof(float) * V_MAX * D));
    }

    inline float * at(int id) {
        return v_b + D * id;
    }

    int add(float * value) {
        int id = size.fetch_add(1);
        for (int i = 0; i < D; i++) {
            v_b[id * D + i] = value[i];
        }
        return id;
    }

    int getSize() {
        return size.load();
    }
};

class ClusterIndex {
public:
    std::vector<std::array<int, LAYER_DEP>> index;

    ClusterIndex(int size) {
        index.resize(size);
    }

    inline int getSize() {
        return index.size();
    }

    inline int at(int i, int lay) {
        return index.at(i).at(lay);
    }
};

class LayIndex {
public:
    std::vector<std::vector<int>> index;
    LayIndex() {
        index.resize(LAYER_DEP, std::vector<int>(V_MAX, 0));
    }

    inline int getLocalId(int layer_id, int global_id) {
        return index.at(layer_id).at(global_id);
    }

    inline void setLocalId(int layer_id, int global_id, int local_id) {
        index.at(layer_id).at(global_id) = local_id;
    }
};

template<int D, int K>
class Layer {
public:
    int layer_id;
    int d_size;
    float* base;
    std::vector<std::vector<int>> neighbors;
    std::vector<std::array<int, K>> nn;
    std::vector<int> next_layer_nearest;
    std::vector<std::atomic<short>>* lock;
    int root;
    bool is_base;

    Layer(int layer_id_, float * base_, int d_size_, std::vector<std::atomic<short>>* lock_, bool is_base_) {

        layer_id = layer_id_;
        d_size = d_size_;
        base = base_;
        lock = lock_;
        is_base = is_base_;
        next_layer_nearest.resize(d_size, -1);
        neighbors.resize(d_size);
        for (int i = 0; i < lock->size(); i++) {
            lock->at(i).store(layer::init);
        }
    }

    inline  bool tryLock(int id) {
        return lock->at(id).compare_exchange_strong(layer::init, layer::layer_thread_id);
    }

    inline void release(int id) {
        lock->at(id).store(layer::init);
    }

    inline float * at(int id) {
        return base + id * D;
    }

    // 要提前clear visited to_visit
    inline int searchOne(float * query, int start, std::pmr::unordered_set<int>& visited, std::vector<int>& to_visit) {
        if(start == -1) {
            start = root;
        }
        int curId = start;
        bool changed = true;
        float cur_dis = vector_dis<D>(query, at(curId));
        while(changed) {
            changed = false;
            auto& neigh = neighbors.at(curId);
            for (int i = 0; i < neigh.size(); i++) {
                if (visited.find(neigh.at(i)) == visited.end()) {
                    to_visit.push_back(neigh.at(i));
                    prefetch<D>(neighbors.at(neigh.at(i)).data());
                }
            }

            for (int i = 0; i < to_visit.size(); i++) {
                visited.insert((to_visit.at(i)));
                float d = vector_dis<D>(query, at(to_visit.at(i)));
                if (d < cur_dis) {
                    cur_dis = d;
                    curId = to_visit.at(i);
                    changed = true;
                }
            }
            to_visit.clear();
        }
        return curId;
    }

    inline searchK(float * query, int * res, int start, int k, std::pmr::unordered_set<int>& visited, std::vector<int>& to_visit) {
        // TODO 要不要挪到solution
        std::priority_queue<std::pair<float, int>> top;
        std::priority_queue<std::pair<float, int>> candidates;

        float cur_dis = vector_dis<D>(at(start), query);
        top.emplace(cur_dis, start);
        candidates.emplace(-cur_dis, start);
        while (candidates.size() > 0) {
            auto &neigh = neighbors.at(candidates.top().second);
            for (int i = 0; i < neigh.size(); i++) {
                if (visited.find(neigh.at(i)) == visited.end()) {
                    to_visit.push_back(neigh.at(i));
                    prefetch<D>(neighbors.at(neigh.at(i)).data());
                }
            }
            candidates.pop();
            if (neigh.size() < 1) continue;

            for (int i = 0; i < to_visit.size(); i++) {
                int id = to_visit.at(i);
                visited.insert(id);
                float d = vector_dis<D>(query, at(id));
                top.emplace(d, id);
                candidates.emplace(-d, id);
            }
            to_visit.clear();

            while (top.size() > k) {
                top.pop();
            }
            if (top.size() == k && top.top().first < -candidates.top().first) {
                for (int rr = 0; rr < k; rr++) {
                    res[rr] = top.top().second;
                    top.pop();
                }
                return;
            }
        }
        for (int rr = 0; rr < k && top.size() > 0; rr++) {
            res[rr] = top.top().second;
            top.pop();
        }
        return;
    }

    void connect(Layer<D,K>* lower, ClusterIndex* cluster_index, LayIndex* lay_index) {
        std::vector<std::pmr::unordered_set<int>> clu(d_size);
        for (int i = 0; i < cluster_index->getSize(); i++) {
            int my_local = lay_index->getLocalId(layer_id, cluster_index->at(i, layer_id));
            int lower_local;
            if (lower->is_base) {
                lower_local = i;
            } else {
                lower_local = lay_index->getLocalId(lower->layer_id, cluster_index->at(i, lower->layer_id));
            }
            clu.at(my_local).insert(lower_local);
        }

        for (int i = 0; i < d_size; i++) {
            int max = 0;
            int max_id = -1;
            for (int low : clu.at(i)) {
                if (lower->neighbors.at(low).size() > max) {
                    max_id = low;
                    max = lower->neighbors.at(low).size();
                }
            }
            next_layer_nearest.at(i) = max_id;
        }
    }

    void setRoot() {
        int max = 0;
        for (int i = 0; i < d_size; i++) {
            if (neighbors.at(i).size() > max) {
                max = neighbors.at(i).size();
                root = i;
            }
        }
    }
};

struct node {
    std::vector<int> neighbors;
    node(int size) {
        neighbors.resize(size);
    }

    node() {
    }

    inline void add_neighbor(int id) {
        neighbors.push_back(id);
    }

    inline int get_neighbor(int id) {
        return neighbors.at(id);
    }

    inline int get_n_size() {
        return neighbors.size();
    }

    inline void clear() {
        neighbors.clear();
    }
};

template<int D, int K>
class Kbrute {
public:
    float * base;
    int d_size;
    int m;
    std::vector<node>* nodes;
    std::vector<int>* id_index;
    std::vector<std::vector<float>>* dis_matric;

    Kbrute(float * base_, int d_size_, int m_, std::vector<node>* nodes_, std::vector<int>* id_index_,
        std::vector<std::vector<float>>* dis_matric_) {
        base = base_;
        d_size = d_size_;
        m = m_;
        nodes = nodes_;
        id_index = id_index_;
        dis_matric = dis_matric_;
    }

    inline float* at(int i) {
        return base + D * (getId(i));
    }

    inline int getId(int i) {
        return id_index->at(i);
    }

    void compute() {
        std::priority_queue<std::pair<float, int>> n_queue;
        int nnn[m];
        for (int i = 0; i < d_size; i++) {
            while (!n_queue.empty()) {
                n_queue.pop();
            }
            for (int j = 0; j < d_size; j++) {
                if (i == j) continue;
                float dis = 0.0f;
                if (i < j) {
                    dis = my_dis(i, j);
                    dis_matric->at(i).at(j) = dis;
                } else {
                    dis = dis_matric->at(j).at(i);
                }
                if (n_queue.size() >= m) {
                    if (dis < n_queue.top().first) {
                        n_queue.emplace(dis, j);
                        n_queue.pop();
                    }
                } else {
                    n_queue.emplace(dis, j);
                }
            }
            for (int x = 0; x < m; x++) {
                nnn[m - 1 - x] = n_queue.top().second;
                n_queue.pop();
            }
            for (int x = 0; x < m; x++) {
                nodes->at(i).add_neighbor(nnn[x]);
            }
        }
    }

    void doMrng() {
        std::vector<int> neigh;
        for (int i = 0; i < d_size; i++) {
            neigh.clear();
            for (int j = 0; j < m; j++) {
                int node_id = nodes->at(i).get_neighbor(j);
                bool add = true;
                float distemp = get_dis(node_id, i);
                for (int k = 0; k < neigh.size(); k++) {
                    if (distemp > get_dis(node_id, neigh.at(k))) {
                        add = false;
                        break;
                    }
                }
                if (add) {
                    neigh.push_back(node_id);
                }
            }
            nodes->at(i).clear();
            for (int zz : neigh) {
                nodes->at(i).add_neighbor(zz);
            }
        }
    }

    inline float my_dis(int a, int b) {
        return vector_dis<D>(at(a), at(b));
    }

    inline float get_dis(int a, int b) {
        if (a < b) {
            return dis_matric->at(a).at(b);
        }
        return dis_matric->at(b).at(a);
    }

    void setKnn(Layer<D, K>* layer_t, LayIndex* layIndex) {
        layer_t->nn.resize(d_size);
        for (int i = 0; i < d_size; i++) {
            for (int j = 0; j < K; j++) {
                layer_t->nn.at(layIndex->getLocalId(layer_t->layer_id, getId(i))).at(j) =
                    layIndex->getLocalId(layer_t->layer_id, getId(nodes->at(i).get_neighbor(j)));
            }
        }
    }

    void addNeighbor(Layer<D, K>* layer_t, LayIndex* layIndex) {
        std::vector<bool> finish(d_size, false);
        std::unordered_set<int> temp_neighbor;
        int count = 0;
        while (count < d_size) {
            for (int i = 0; i < d_size; i++) {
                if (finish.at(i)) continue;
                if (!layer_t->tryLock(layIndex->getLocalId(layer_t->layer_id, getId(i)))) continue;

                int layer_local_i = layIndex->getLocalId(layer_t->layer_id, getId(i));
                temp_neighbor.insert(layer_t->neighbors.at(layer_local_i).begin(), layer_t->neighbors.at(layer_local_i).end());
                for (int j = 0; j < nodes->at(i).get_n_size(); j++) {
                    int layer_local_j = layIndex->getLocalId(layer_t->layer_id, getId(nodes->at(i).get_neighbor(j)));
                    temp_neighbor.insert(layer_local_j);
                }
                layer_t->neighbors.at(layer_local_i).assign(temp_neighbor.begin(), temp_neighbor.end());

                count++;
                temp_neighbor.clear();
                finish.at(i) = true;
                layer_t->release(layIndex->getLocalId(layer_t->layer_id, getId(i)));
            }
        }
    }

    void addBaseNeighbor(Layer<D, K>* layer_t) {
        std::vector<bool> finish(d_size, false);
        std::unordered_set<int> temp_neighbor;
        int count = 0;
        while (count < d_size) {
            for (int i = 0; i < d_size; i++) {
                if (finish.at(i)) continue;
                if (!layer_t->tryLock(getId(i)))continue;

                int layer_local_i = getId(i);
                temp_neighbor.insert(layer_t->neighbors.at(layer_local_i).begin(), layer_t->neighbors.at(layer_local_i).end());
                for (int j = 0; j < nodes->at(i).get_n_size(); j++) {
                    int layer_local_j = getId(nodes->at(i).get_neighbor(j));
                    temp_neighbor.insert(layer_local_j);
                }
                layer_t->neighbors.at(layer_local_i).assign(temp_neighbor.begin(), temp_neighbor.end());

                count++;
                temp_neighbor.clear();
                finish.at(i) = true;
                layer_t->release(layer_local_i);
            }
        }
    }
};

template<int D, int K>
void build(Layer<D, K>* layer_t, ClusterIndex* clusterIndex, LayIndex* layIndex, VBase<D>* vb) {
    int out = std::min(layer_t->d_size - 1, MAX_M);
    std::vector<node> nodes(layer_t->d_size, out);
    std::vector<std::vector<float>> dis_matric(layer_t->d_size);
    for (int i = 0; i < layer_t->d_size; i++) {
        dis_matric.at(i).resize(layer_t->d_size);
    }
    std::vector<int> cluster_ids;
    std::pmr::unordered_set<int> ids;
    for (int j = 0; j < clusterIndex->getSize(); j++) {
        int this_global_id = clusterIndex->at(j, layer_t->layer_id);
        ids.insert(this_global_id);
    }
    cluster_ids.assign(ids.begin(), ids.end());
    Kbrute<D, K> kb(vb->at(0), layer_t->d_size, out, &nodes, &cluster_ids, &dis_matric);
    kb.compute();
    kb.setKnn(layer_t, layIndex);
    kb.doMrng();
    kb.addNeighbor(layer_t, layIndex);
    layer_t->setRoot();
}

template<int D, int K>
void crossBuild(Layer<D, K>* layer_t, ClusterIndex* clusterIndex, LayIndex* layIndex, VBase<D>* vb, float* real_base,
    Layer<D, K>* layer_cross) {
    std::vector<node> nodes;
    std::vector<std::vector<float>> dis_matric_;
    std::vector<int> cluster_ids;
    std::unordered_set<int> cluster_ids_set;
    std::priority_queue<std::pair<float, int>> cross;

    if (layer_t->is_base) {
        for (int i = 0; i < layer_cross->d_size; i++) {
            cluster_ids_set.clear();
            for (int j = 0; j < clusterIndex->getSize(); j++) {
                // i 是 cross_local_id
                int cross_global_id = clusterIndex->at(j, layer_cross->layer_id);
                int this_global_id = j;
                if (layIndex->getLocalId(layer_cross->layer_id, cross_global_id) == i) {
                    cluster_ids_set.insert(this_global_id);
                }
            }
            // select cross
            for (int near : layer_cross->nn.at(i)) {
                for (int j = 0; j < clusterIndex->getSize(); j++) {
                    // near 是cross_local_id
                    int cross_global_id = clusterIndex->at(j, layer_cross->layer_id);
                    int this_global_id = j;
                    if (layIndex->getLocalId(layer_cross->layer_id, cross_global_id) == near) {
                        float dis = vector_dis<D>(layer_cross->at(i), layer_t->at(this_global_id));
                        if (cross.size() >= OVER_LAP && dis > cross.top().first) continue;
                        cross.emplace(dis, this_global_id);
                        if (cross.size() > OVER_LAP) {
                            cross.pop();
                        }
                    }
                }
                while (!cross.empty()) {
                    cluster_ids_set.insert(cross.top().second);
                    cross.pop();
                }
            }

            cluster_ids.assign(cluster_ids_set.begin(), cluster_ids_set.end());
            int k_s = std::min<int>(cluster_ids.size() - 1, MAX_M);
            nodes.clear();
            nodes.resize(cluster_ids.size(), k_s);
            dis_matric_.resize(cluster_ids.size());
            for (auto& ddd : dis_matric_) {
                ddd.resize(cluster_ids.size());
            }

            Kbrute<D, K> kb(real_base, cluster_ids.size(), k_s, &nodes, &cluster_ids, &dis_matric_);
            kb.compute();
            kb.doMrng();
            kb.addBaseNeighbor(layer_t);
        }
    } else {
        for (int i = 0; i < layer_cross->d_size; i++) {
            cluster_ids_set.clear();
            for (int j = 0; j < clusterIndex->getSize(); j++) {
                // i 是 cross_local_id
                int cross_global_id = clusterIndex->at(j, layer_cross->layer_id);
                int this_global_id = clusterIndex->at(j, layer_t->layer_id);
                if (layIndex->getLocalId(layer_cross->layer_id, cross_global_id) == i) {
                    cluster_ids_set.insert(this_global_id);
                }
            }
            // select cross
            for (int near : layer_cross->nn.at(i)) {
                for (int j = 0; j < clusterIndex->getSize(); j++) {
                    // near 是cross_local_id
                    int cross_global_id = clusterIndex->at(j, layer_cross->layer_id);
                    int this_global_id = clusterIndex->at(j, layer_t->layer_id);

                    if (layIndex->getLocalId(layer_cross->layer_id, cross_global_id) == near) {
                        float dis = vector_dis<D>(layer_cross->at(i), layer_t->at(this_global_id));
                        if (cross.size() >= OVER_LAP && dis > cross.top().first) continue;
                        cross.emplace(dis, this_global_id);
                        if (cross.size() > OVER_LAP) {
                            cross.pop();
                        }
                    }
                }
                while (!cross.empty()) {
                    cluster_ids_set.insert(cross.top().second);
                    cross.pop();
                }
            }

            cluster_ids.assign(cluster_ids_set.begin(), cluster_ids_set.end());
            int k_s = std::min<int>(cluster_ids.size() - 1, MAX_M);
            nodes.clear();
            nodes.resize(cluster_ids.size(), k_s);
            dis_matric_.resize(cluster_ids.size());
            for (auto& ddd : dis_matric_) {
                ddd.resize(cluster_ids.size());
            }

            Kbrute<D, K> kb(vb->at(0), cluster_ids.size(), k_s, &nodes, &cluster_ids, &dis_matric_);
            kb.compute();
            kb.doMrng();
            kb.addBaseNeighbor(layer_t);
        }
    }
}

template<int D, int K>
Layer<D, K>* create(int id, VBase<D>* vb, LayIndex* layIndex, std::pmr::unordered_set<int>& ids) {
    int local_id = 0;
    float* lb = static_cast<float*>(aligned_alloc(CACHE_LINE, sizeof(float) * D * ids.size()));
    for (int n : ids) {
        layIndex->setLocalId(id, n, local_id);
        std::memcpy(lb + D * local_id, vb->at(n), D * sizeof(float));
        local_id++;
    }
    std::vector<std::atomic<short>>* lock = new std::vector<std::atomic<short>>[ids.size()];
    Layer<D, K>* layer = new Layer<D, K>(id, lb, ids.size(), lock, false);
    return layer;
}

template<int D, int K>
std::vector<Layer<D, K>*> createAll(VBase<D>* vb, float * real_base, ClusterIndex* clusterIndex) {
    std::vector<Layer<D, K>*> res;

    std::unordered_set<int> ids0(64);
    std::unordered_set<int> ids1(2048);
    std::unordered_set<int> ids2(32 * 1024);
    for (int i = 0; i < clusterIndex->getSize(); i++) {
        ids0.insert(clusterIndex->at(i, 0));
        ids1.insert(clusterIndex->at(i, 1));
        ids2.insert(clusterIndex->at(i, 2));
    }

    LayIndex layIndex;
    Layer<<D,K>* layer0 = create<D, K>(0, vb, &layIndex, ids0);
    build(layer0, clusterIndex, &layIndex, vb);
    layer0.setRoot();

    Layer<<D,K>* layer1 = create<D, K>(1, vb, &layIndex, ids1);
    build(layer1, clusterIndex, &layIndex, vb);
    layer0.connect(layer1, clusterIndex, &layIndex);

    Layer<<D,K>* layer2 = create<D, K>(2, vb, &layIndex, ids2);
    crossBuild(layer2, clusterIndex, &layIndex, vb, real_base, layer0);
    layer1.connect(layer2, clusterIndex, &layIndex);

    std::vector<std::atomic<short>>* lock = new std::vector<std::atomic<short>>(clusterIndex->getSize());
    Layer<<D,K>* layer_base = new Layer<D, K>(3, real_base, clusterIndex->getSize(), lock, true);
    crossBuild(layer_base, clusterIndex, &layIndex, vb, real_base, layer1);
    layer2.connect(layer_base, clusterIndex, &layIndex);

    res.push_back(layer0);
    res.push_back(layer1);
    res.push_back(layer2);
    res.push_back(layer_base);
    return res;
}


#endif //LAYER_H
