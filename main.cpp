#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

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
std::vector<std::vector<float>> readFvecs(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<std::vector<float>> vectors;
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
            vectors.push_back(vector);
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
    try {
        {
            std::string filename = "/root/dataset/newnew/gist/gist_base.fvecs";
            std::vector<std::vector<float>> vectors = readFvecs(filename);

            std::cout << "Read " << vectors.size() << " vectors with dimension " << vectors[0].size() << std::endl;

            // 打印前几个向量
            for (size_t i = 0; i < std::min<size_t>(5, vectors.size()); ++i) {
                std::cout << "Vector " << i << ": ";
                for (float value : vectors[i]) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
            }
        }

        {
            std::string filename = "/root/dataset/newnew/gist/gist_groundtruth.ivecs";
            std::vector<std::vector<int>> vectors = readIvecs(filename);

            std::cout << "Read " << vectors.size() << " vectors with dimension " << vectors[0].size() << std::endl;

            // 打印前几个向量
            for (size_t i = 0; i < std::min<size_t>(5, vectors.size()); ++i) {
                std::cout << "Vector " << i << ": ";
                for (int value : vectors[i]) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
            }

        }




    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
//
// int main() {
//     try {
//         std::vector<std::vector<float>> vectors = {
//             {1.0, 2.0, 3.0},
//             {4.0, 5.0, 6.0},
//             {7.0, 8.0, 9.0}
//         };
//
//         std::string filename = "output.fvecs";
//         writeFvecs(filename, vectors);
//
//         std::cout << "Wrote " << vectors.size() << " vectors with dimension " << vectors[0].size() << " to " << filename << std::endl;
//     } catch (const std::exception& e) {
//         std::cerr << "Exception: " << e.what() << std::endl;
//         return 1;
//     }
//
//     return 0;
// }