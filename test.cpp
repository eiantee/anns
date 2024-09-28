//
// Created by caiqibang 00906710 on 2024/9/28.
//

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <iostream>
#include <vector>
#include <faiss/IndexFlat.h>
#include <faiss/Clustering.h>
#include <faiss/IndexHNSW.h>
#include <faiss/MetricType.h>
#include <faiss/impl/FaissAssert.h>
#include <faiss/impl/AuxIndexStructures.h>
#include <faiss/VectorTransform.h>
#include <faiss/Clustering.h>

#include "hnsw.h"

// 写入 fvecs 文件
void writeFvecs(const std::string& filename, const std::vector<std::vector<float>>& vectors) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    if (vectors.empty()) {
        throw std::runtime_error("No vectors to write");
    }

    uint32_t dimension = static_cast<uint32_t>(vectors[0].size());

    // 写入文件头，维度
    file.write(reinterpret_cast<const char*>(&dimension), sizeof(dimension));
    if (!file) {
        throw std::runtime_error("Could not write dimension to file: " + filename);
    }

    // 写入向量数据
    for (const auto& vector : vectors) {
        if (vector.size() != dimension) {
            throw std::runtime_error("Inconsistent vector dimensions");
        }
        file.write(reinterpret_cast<const char*>(vector.data()), dimension * sizeof(float));
        if (!file) {
            throw std::runtime_error("Error writing to file: " + filename);
        }
    }
}

// 读取 fvecs 文件
std::vector<float> readFvecs(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<float> vectors;
    uint32_t dimension;

    // 读取文件头，获取维度
    file.read(reinterpret_cast<char*>(&dimension), sizeof(dimension));
    if (!file) {
        throw std::runtime_error("Could not read dimension from file: " + filename);
    }

    // 读取每个向量
    while (file) {
        std::vector<float> vector(dimension);
        file.read(reinterpret_cast<char*>(vector.data()), dimension * sizeof(float));
        if (file) {
            vectors.insert(vectors.end(), vector.begin(), vector.end());
        }
    }

    if (file.bad()) {
        throw std::runtime_error("Error reading file: " + filename);
    }

    return vectors;
}


// 读取 Ivecs 文件
std::vector<std::vector<int>> readIvecs(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<std::vector<int>> vectors;
    uint32_t dimension;

    // 读取文件头，获取维度
    file.read(reinterpret_cast<char*>(&dimension), sizeof(dimension));
    if (!file) {
        throw std::runtime_error("Could not read dimension from file: " + filename);
    }

    // 读取每个向量
    while (file) {
        std::vector<int> vector(dimension);
        file.read(reinterpret_cast<char*>(vector.data()), dimension * sizeof(float));
        if (file) {
            vectors.push_back(vector);
        }
    }

    if (file.bad()) {
        throw std::runtime_error("Error reading file: " + filename);
    }

    return vectors;
}

int main() {

    std::string filename = "/root/dataset/newnew/gist/gist_base.fvecs";
    std::vector<float> vectors = readFvecs(filename);

    // 数据和参数初始化
    constexpr int d = 960; // 维度
    int nb = vectors.size() /960 ; // 数据点数量
    int ncentroids = 16; // 聚类中心的数量



    Hnsw<d>* hsh = new Hnsw<d>(vectors.data(), vectors.size(), 20, 40, 200, 200);
    hsh->build();


    // 定义索引和聚类器
    faiss::IndexFlatL2 index(d);  // 使用L2距离
    faiss::Clustering clustering(d, ncentroids);

    // 配置聚类参数
    faiss::ClusteringParameters params;
    params.niter = 20;
    params.verbose = true;

    // faiss::IndexHNSWFlat indxxex(d, 32);
    // indxxex.add(nb, vectors.data());

    // 执行聚类
    clustering.train(nb, vectors.data(), index);

    // 获取聚类结果
    std::vector<float> centroids(ncentroids * d);
    index.reconstruct_n(0, ncentroids, centroids.data());

    return 0;

}