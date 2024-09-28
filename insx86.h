//
// Created by caiqibang 00906710 on 2024/9/26.
//

#ifndef INSX86_H
#define INSX86_H
#include <immintrin.h>

static inline void prefetch1h(const void *data) {
    _mm_prefetch(data, _MM_HINT_T0);
}

static inline void prefetch3h(const void *data) {
    _mm_prefetch(data, _MM_HINT_T2);
}

static inline void prefetch1l(const void *data) {
    _mm_prefetch(data, _MM_HINT_T0);
}

static inline void prefetch3l(const void *data) {
    _mm_prefetch(data, _MM_HINT_T2);
}

template<int S>
inline void sub(const float* a, const float* b, float* result) {
    sub<4>(a, b, result);
    sub<S - 4>(a + 4, b + 4, result + 4);
}

template<>
inline void sub<4>(const float* a, const float* b, float* result) {
    __m128 va = _mm_loadu_ps(a);
    __m128 vb = _mm_loadu_ps(b);
    va = _mm_sub_ps(va, vb);
    _mm_storeu_ps(result, va);
}

template<>
inline void sub<3>(const float* a, const float* b, float* result) {
    for (int i = 0; i < 3; i++) {
        result[i] = a[i] - b[i];
    }
}

template<>
inline void sub<2>(const float* a, const float* b, float* result) {
    for (int i = 0; i < 2; i++) {
        result[i] = a[i] - b[i];
    }
}

template<>
inline void sub<1>(const float* a, const float* b, float* result) {
    result[0] = a[0] - b[0];
}

template<int S>
inline __m128 l2f(const float* a, const float* b) {
    __m128 va = l2f<4>(a, b);
    __m128 vb = l2f<S - 4>(a + 4, b + 4);
    return _mm_add_ps(va, vb);
}

template<>
inline __m128 l2f<4>(const float* a, const float* b) {
    __m128 va = _mm_loadu_ps(a);
    __m128 vb = _mm_loadu_ps(b);
    va = _mm_sub_ps(va, vb);
    return _mm_mul_ps(va, va);
}

template<>
inline __m128 l2f<3>(const float* a, const float* b) {
    __m128 va = _mm_set_ps(a[0], a[1], a[2], 0.0f);
    __m128 vb = _mm_set_ps(b[0], b[1], b[2], 0.0f);
    va = _mm_sub_ps(va, vb);
    return _mm_mul_ps(va, va);
}

template<>
inline __m128 l2f<2>(const float* a, const float* b) {
    __m128 va = _mm_set_ps(a[0], a[1], 0.0f, 0.0f);
    __m128 vb = _mm_set_ps(b[0], b[1], 0.0f, 0.0f);
    va = _mm_sub_ps(va, vb);
    return _mm_mul_ps(va, va);
}

template<>
inline __m128 l2f<1>(const float* a, const float* b) {
    __m128 va = _mm_set_ps(a[0], 0.0f, 0.0f, 0.0f);
    __m128 vb = _mm_set_ps(b[0], 0.0f, 0.0f, 0.0f);
    va = _mm_sub_ps(va, vb);
    return _mm_mul_ps(va, va);
}

template<int D>
float vector_dis(const float* a, const float* b) {
    __m128 vec = l2f<D>(a, b);
    // 第一步：将高 64 位和低 64 位的元素相加
    __m128 sum1 = _mm_add_ps(vec, _mm_movehl_ps(vec, vec));

    // 第二步：将 128 位向量的低 32 位和高 32 位的元素相加
    __m128 sum2 = _mm_add_ss(sum1, _mm_shuffle_ps(sum1, sum1, 1));

    // 将结果转换为 float
    return _mm_cvtss_f32(sum2);
}


template<typename T>
void assert_equal_all(T * a, T * b, int size) {
    for (int i = 0; i < size; i++) {
        assert(abs(a[i] - b[i]) < 0.000001f);
    }
}

template<typename T>
void assert_equal(T a, T b) {
    assert(abs(a - b) < 0.000001f);
}

void assert_arch() {
    float a[1027];
    float b[1027];
    for (int i = 0; i < 1027; i++) {
        a[i] = i % 11;
        b[i] = i % 17;
    }

    float l2 = 0.0f;
    for (int i = 0; i < 1027; i++) {
        l2 += (a[i] - b[i]) * (a[i] - b[i]);
    }
    assert_equal(l2, vector_dis<1027>(a, b));

    float c[1027];
    sub<1027>(a, b, c);
    float d[1027];
    for (int i = 0; i < 1027; i++) {
        d[i] = a[i] - b[i];
    }
    assert_equal_all(c, d, 1027);
}

#endif //INSX86_H
