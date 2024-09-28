//
// Created by caiqibang 00906710 on 2024/9/28.
//

#ifndef HNSW_H
#define HNSW_H

#include <atomic>
#include <mutex>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "config.h"

#define MAX_LAYER 4
#define L_INIT -1
#define NEI_LENGTH 256

class Neighbor {
public:
    // std::atomic<short> lock;
    int size;
    int neighbor_ids[(NEI_LENGTH - sizeof(int)) / sizeof(int)];

    inline short getSize() {
        return size;
    }

    inline int getNeighborId(int index) {
        return neighbor_ids[index];
    }

    inline void add(int id) {
        neighbor_ids[size] = id;
        size++;
    }

    inline void clear() {
        size = 0;
    }
};

struct alignas(CACHE_LINE) LockWrapper {
    std::atomic<int> lock{L_INIT};

    void release(int thread_id) {
        int expect = thread_id;
        lock.compare_exchange_strong(expect, L_INIT);
    }

    bool tryLock(int thread_id) {
        int desire = L_INIT;
        return lock.compare_exchange_strong( desire, thread_id);
    }


};

class NeighborHandler {
public:
    Neighbor neighbors[MAX_LAYER];
    inline short getSize(int layer) {
        return neighbors[layer].getSize();
    }

    inline int getNeighborId(int index, int layer) {
        return neighbors[layer].getNeighborId(index);
    }

    inline void add(int id, int layer) {
        neighbors[layer].add(id);
    }

    inline void clear(int layer) {
        neighbors[layer].clear();
    }
};

template<int D>
class Hnsw {
public:
    static Hnsw<D>* global_hnsw = nullptr;
    const float* base;
    int size;
    int max_out{15};
    int max_out_0{31};
    int ef_add;
    int ef_search;
    int root;

    bool * added;
    LockWrapper* locks;
    int * labels;

    // 为了debug 展示 可以临时改为vector
    NeighborHandler* neighbors;


    Hnsw(const float* base, int size, int max_out, int max_out_0, int ef_add, int ef_search) :
    base(base), size(size), max_out(max_out), max_out_0(max_out_0), ef_add(ef_add), ef_search(ef_search) {
        void * address = aligned_alloc(CACHE_LINE, sizeof(NeighborHandler) * size);
        neighbors = new (address) NeighborHandler[size];


        locks = new LockWrapper[size];
        labels = new int[size]{-1};

        added = new bool[size]{false};
    }

    static void * run_add_base(void * args) {
        int thread_id = *((int *) args);
        my_thread_id = thread_id;
        setAffinity(thread_id);

        // thread segment
        int segment = (Hnsw<D>::global_hnsw->size - 1) / THREADS + 1;
        int start_index = segment * my_thread_id;
        int end_index = std::min<int>(start_index + segment, Hnsw<D>::global_hnsw->size);
        int left = end_index - start_index;

        // 获取锁失败的时候 要记住上回 需要锁哪个， 然后提前锁
        std::unordered_map<int, std::unordered_set<int>> lock_map;
        while (true) {
            for (int i = start_index; i < end_index; i++) {
                if (Hnsw<D>::global_hnsw->added[i]) continue;
                if (Hnsw<D>::global_hnsw->add_base(i, lock_map) == -1) continue;
                Hnsw<D>::global_hnsw->added[i] = true;
                left--;
                if (left < 1) {
                    return nullptr;
                }
            }
        }
    }

    ~Hnsw() {
        free(neighbors);
    }

    inline const float* at(int id) {
        return base + D * id;
    }

    inline Neighbor get_neighbor(int id, int level) {
        return neighbors[id].neighbors[level];
    }

    void build() {
        // 随机
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, size - 1);

        root = dis(gen);
        added[root] = true;

        // root 层 32个 MAX_LAYER - 1
        std::unordered_set<int> node_l1;
        while(node_l1.size() < 31) {
            int id = dis(gen);
            if (id == root) continue;
            node_l1.insert(id);

            added[id] = true;
        }
        for (int node_id : node_l1) {
            add(node_id, MAX_LAYER - 1);
        }

        // 二 层 1024个 MAX_LAYER - 2
        std::unordered_set<int> node_l2;
        while(node_l2.size() < 1024 - 32) {
            int id = dis(gen);
            if (id == root) continue;
            if (node_l1.find(id) != node_l1.end()) continue;
            node_l2.insert(id);

            added[id] = true;
        }
        for (int node_id : node_l2) {
            add(node_id, MAX_LAYER - 2);
        }

        // 三 层 10000 MAX_LAYER - 3
        std::unordered_set<int> node_l3;
        while(node_l3.size() < 10000 - 1024) {
            int id = dis(gen);
            if (id == root) continue;
            if (node_l1.find(id) != node_l1.end()) continue;
            if (node_l2.find(id) != node_l2.end()) continue;
            node_l3.insert(id);

            added[id] = true;
        }
        for (int node_id : node_l3) {
            add(node_id, MAX_LAYER - 3);
        }

        global_hnsw = this;

        // 多线程 add_base
        pthread_t threads[THREADS - 1];
        int* thread_ids[THREADS - 1];
        for(int i = 0; i < THREADS - 1; i++) {
            thread_ids[i] = (int*) malloc(sizeof(int));
            *thread_ids[i] = i + 1;

            pthread_create(&threads[i], NULL, run_add_base, (void*) thread_ids[i]);
        }
    }

    void add(int node_id, const int add_level) {
        // 通过cluster center 一开始就创建最大的 level
        int level_enter = root;

        if (add_level < MAX_LAYER - 1) {
            float curdist = vector_dis<D>(at(node_id), at(root));
            for (int level = MAX_LAYER - 1; level > add_level; level--) {
                bool changed = true;
                while (changed) {
                    changed = false;
                    // TODO lock
                    Neighbor neighbor = get_neighbor(level_enter, level);
                    for (int i = 0; i < neighbor.getSize(); i++) {
                        float dis = vector_dis<D>(at(node_id), at(neighbor.getNeighborId(i)));
                        if (dis < curdist) {
                            curdist = dis;
                            level_enter = neighbor.getNeighborId(i);
                            changed = true;
                        }
                    }
                }
            }
        }
        std::unordered_set<int> visited;
        std::vector<int> to_visit;
        std::priority_queue<std::pair<float, int>> top_can;
        for (int level = add_level; level >= 0; level--) {
            visited.clear();
            to_visit.clear();
            while (!top_can.empty()) {
                top_can.pop();
            }
            search_dis(at(node_id), top_can, level_enter, level, ef_add, visited, to_visit);
            top_can.emplace(vector_dis<D>(at(node_id), at(level_enter)), level_enter);
            level_enter = mutuallyConnect(node_id, top_can, level);
        }
    }

    int mutuallyConnect(int node_id, std::priority_queue<std::pair<float, int>>& top_can, int level) {
        int Mcurmax = (level == 0 ? max_out_0 : max_out);
        prune(top_can, Mcurmax);

        std::vector<int> selectedNeighbors;
        selectedNeighbors.reserve(Mcurmax);
        while (!top_can.empty()) {
            selectedNeighbors.push_back(top_can.top().second);
            top_can.pop();
        }

        int next_closest_entry_point = selectedNeighbors.back();

        // TODO 互相更新怎么锁
        std::priority_queue<std::pair<float, int>> nei_top_can;
        for (int idx : selectedNeighbors) {
            Neighbor nei_neighbor = get_neighbor(idx, level);
            if (nei_neighbor.getSize() < Mcurmax) {
                nei_neighbor.add(node_id);
            } else {
                float dis = vector_dis<D>(at(node_id), at(idx));
                nei_top_can.emplace(dis, node_id);
                for (int nn = 0; nn < nei_neighbor.getSize(); nn++) {
                    float ndis = vector_dis<D>(at(nn), at(idx));
                    nei_top_can.emplace(ndis, nn);
                }
                prune(nei_top_can, Mcurmax);
                nei_neighbor.clear();
                // 不一定要互相连
                while (!nei_top_can.empty()) {
                    nei_neighbor.add(nei_top_can.top().second);
                    nei_top_can.pop();
                }
            }
        }
        return next_closest_entry_point;
    }

    void prune(std::priority_queue<std::pair<float, int>>& top_can, int max) {
        if (top_can.size() < max) {
            return;
        }

        std::priority_queue<std::pair<float, int>> queue_closest;
        std::vector<std::pair<float, int>> return_list;
        while (!top_can.empty()) {
            queue_closest.emplace(-top_can.top().first, top_can.top().second);
            top_can.pop();
        }

        while (!queue_closest.empty()) {
            if (return_list.size() >= max)
                break;
            std::pair<float, int> current_pair = queue_closest.top();
            float dist_to_query = -current_pair.first;
            queue_closest.pop();
            bool good = true;

            for (std::pair<float, int> second_pair : return_list) {
                float curdist = vector_dis<D>(at(current_pair.second), at(second_pair.first));
                if (curdist < dist_to_query) {
                    good = false;
                    break;
                }
            }
            if (good) {
                return_list.push_back(current_pair);
            }
        }

        for (auto current_pair : return_list) {
            top_can.emplace(-current_pair.first, current_pair.second);
        }
    }

    inline void search_dis(const float * query, std::priority_queue<std::pair<float, int>>& top_can, int start, int level, int ef,
        std::unordered_set<int>& visited, std::vector<int>& to_visit) {
        std::priority_queue<std::pair<float, int>> candidates;

        float cur_dis = vector_dis<D>(at(start), query);
        top_can.emplace(cur_dis, start);
        candidates.emplace(-cur_dis, start);
        visited.insert(start);
        while (!candidates.empty()) {
            auto neigh = get_neighbor(candidates.top().second, level);
            candidates.pop();
            for (int i = 0; i < neigh.getSize(); i++) {
                if (visited.find(neigh.getNeighborId(i)) == visited.end()) {
                    to_visit.push_back(neigh.getNeighborId(i));
                    prefetchLocalL(at(neigh.getNeighborId(i)), D * sizeof(float));
                }
            }
            if (neigh.getSize() < 1) continue;
            for (int i = 0; i < to_visit.size(); i++) {
                int id = to_visit.at(i);
                visited.insert(id);
                float d = vector_dis<D>(query, at(id));
                top_can.emplace(d, id);
                candidates.emplace(-d, id);
            }
            to_visit.clear();

            if (candidates.empty()) return;
            while (top_can.size() > ef) {
                top_can.pop();
            }
            if (top_can.size() == ef && top_can.top().first < -candidates.top().first) {
                return;
            }
        }
    }

    int add_base(int node_id, std::unordered_map<int, std::unordered_set<int>>& label_map) {
        // 通过cluster center 一开始就创建最大的 level
        int level_enter = root;
        int this_lock_label = -1;

        auto label_set = label_map.find(node_id);

        if (label_set != label_map.end()) {
            // 提前锁
            for (int i : label_map[node_id]) {
                if (!locks[i].tryLock(my_thread_id)) {
                    for (int j : label_map[node_id]) {
                        locks[i].release(my_thread_id);
                    }
                    return -1;
                }
            }
        }

        float curdist = vector_dis<D>(at(node_id), at(root));
        for (int level = MAX_LAYER - 1; level > 0; level--) {
            bool changed = true;
            while (changed) {
                changed = false;
                Neighbor neighbor = get_neighbor(level_enter, level);
                for (int i = 0; i < neighbor.getSize(); i++) {
                    float dis = vector_dis<D>(at(node_id), at(neighbor.getNeighborId(i)));
                    if (dis < curdist) {
                        curdist = dis;
                        level_enter = neighbor.getNeighborId(i);
                        changed = true;
                    }
                }
            }

            if (level == 2) {
                labels[node_id] = level_enter;
            }
        }
        // TODO lock 根据label 加锁， label  就是 2 层节点的id, 从一层节点进来就添加 // 1024 个锁


        std::unordered_set<int> visited;
        std::vector<int> to_visit;
        std::priority_queue<std::pair<float, int>> top_can;

        search_dis(at(node_id), top_can, level_enter, 0, ef_add, visited, to_visit);
        top_can.emplace(vector_dis<D>(at(node_id), at(level_enter)), level_enter);

        // 根据 visited 加锁 锁全部

        std::unordered_set<int> label_to_lock;
        for (int l : visited) {
            label_to_lock.insert(labels[l]);
        }

        if (label_set != label_map.end()) {
            for (int i : label_set->second) {
                label_to_lock.erase(i);
            }
        }

        for (int la : label_to_lock) {
            label_set->second.insert(la);
        }

        for (int la : label_to_lock) {
            if (!locks[la].tryLock(my_thread_id)) {
                for (int xx : label_to_lock) {
                    locks[xx].release(my_thread_id);
                }
                return -1;
            }
        }

        level_enter = mutuallyConnect(node_id, top_can, 0);

        for (int i : label_set->second) {
            locks[i].release(my_thread_id);
        }


        return 0;
    }

    inline void search_base(float * query, int * res, int start, int k, float cur_dis) {

        // TODO 要不要挪到solution
        std::unordered_set<int> visited;
        std::vector<int> to_visit;
        std::priority_queue<std::pair<float, int>> top;
        std::priority_queue<std::pair<float, int>> candidates;

        top.emplace(cur_dis, start);
        candidates.emplace(-cur_dis, start);
        visited.insert(start);
        while (candidates.size() > 0) {
            auto &neigh = get_neighbor(candidates.top().second, 0);
            for (int i = 0; i < neigh.getSize(); i++) {
                if (visited.find(neigh.getNeighborId(i)) == visited.end()) {
                    to_visit.push_back(neigh.getNeighborId(i));
                    prefetchLocalL(at(neigh.getNeighborId(i)), D * sizeof(float));
                }
            }
            candidates.pop();
            if (neigh.getSize() < 1) continue;

            for (int i = 0; i < to_visit.size(); i++) {
                int id = to_visit.at(i);
                visited.insert(id);
                float d = vector_dis<D>(query, at(id));
                top.emplace(d, id);
                candidates.emplace(-d, id);
            }
            to_visit.clear();

            while (top.size() > ef_search) {
                top.pop();
            }
            if (top.size() == ef_search && top.top().first < -candidates.top().first) {
                while (top.size() > k) {
                    top.pop();
                }
                for (int rr = 0; rr < k; rr++) {
                    res[k - rr - 1] = top.top().second;
                    top.pop();
                }
                return;
            }
        }
        while (top.size() > k) {
            top.pop();
        }
        for (int rr = std::min<int>(k, top.size()) - 1; rr >= 0; rr--) {
            res[rr] = top.top().second;
            top.pop();
        }
    }

    inline void search(float * query, int * res, int k) {
        int level_enter = root;
        float curdist = vector_dis<D>(query, at(root));
        for (int level = MAX_LAYER - 1; level > 0; level--) {
            bool changed = true;
            while (changed) {
                changed = false;
                Neighbor* neighbor = get_neighbor(level_enter, level);
                for (int i = 0; i < neighbor->getSize(); i++) {
                    float dis = vector_dis<D>(query, at(neighbor->getNeighborId(i)));
                    if (dis < curdist) {
                        curdist = dis;
                        level_enter = neighbor->getNeighborId(i);
                        changed = true;
                    }
                }
            }
        }
        search_base(query, res, level_enter, k, curdist);
    }
};

#endif //HNSW_H
