//
// Created by caiqibang 00906710 on 2024/9/27.
//

#ifndef CPP_SOLUTION_H
#define CPP_SOLUTION_H
#include <vector>

class Solution {
public:
    void build(int d, const std::vector<float>& base);
    void search(const std::vector<float>& query, int* res);
};

#endif //SOLUTION_H
