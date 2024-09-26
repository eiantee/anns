//
// Created by caiqibang 00906710 on 2024/9/26.
//

#ifndef INSX86_H
#define INSX86_H
#include <immintrin.h>

template<int D>
inline void addf(double *a, float *b) {
    for (int i = 0; i < D; i++) {
        a[i] += b[i];
    }
}

template<int D>
inline void scaled(double * a, double b) {
    for (int i = 0; i < D; i++) {
        a[i] *= b;
    }
}

template<int D>
inline float l2f(const float* a, const float* b) {
    float ret = 0.0f;
    for (int i = 0; i < D; i++) {
        float t = a[i] - b[i];
        ret += t * t;
    }
    return ret;
}

template<int D>
float vector_dis(const float* a, const float* b) {
    return l2f<D>(a, b);
}

#endif //INSX86_H
