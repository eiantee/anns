//
// Created by caiqibang 00906710 on 2024/9/26.
//

#ifndef INSARM_H
#define INSARM_H
#include <arm_neon.h>
#include <cassert>
#include <cmath>

static inline void prefetch1h(const void *data) {
    __asm__ __volatile__(
        "prfm PLDL1KEEP, [%[data]]    \n\t"
        :: [data] "r" (data));
}

static inline void prefetch3h(const void *data) {
    __asm__ __volatile__(
        "prfm PLDL3KEEP, [%[data]]    \n\t"
        :: [data] "r" (data));
}

static inline void prefetch1l(const void *data) {
    __asm__ __volatile__(
        "prfm PLDL1STRM, [%[data]]    \n\t"
        :: [data] "r" (data));
}

static inline void prefetch3l(const void *data) {
    __asm__ __volatile__(
        "prfm PLDL3STRM, [%[data]]    \n\t"
        :: [data] "r" (data));
}

//https://developer.arm.com/documentation/den0018/a/NEON-Intrinsics-Reference/Arithmetic/VSUB?lang=en
template<int S>
inline void sub(const float* a, const float* b, float* result) {
    sub<4>(a, b, result);
    sub<S - 4>(a + 4, b + 4, result + 4);
}

template<>
inline void sub<4>(const float* a, const float* b, float* result) {
    float32x4_t va = vld1q_f32(a);
    float32x4_t vb = vld1q_f32(b);
    va = vsubq_f32(va, vb);
    vst1q_f32(result, va);
}

template<>
inline void sub<3>(const float* a, const float* b, float* result) {
    float32x4_t va = {a[0], a[1], a[2], 0.0f};
    float32x4_t vb = {b[0], b[1], b[2], 0.0f};
    va = vsubq_f32(va, vb);
    vst1q_lane_f32(result, va, 0);
    vst1q_lane_f32(result + 1, va, 1);
    vst1q_lane_f32(result + 2, va, 2);
}

template<>
inline void sub<2>(const float* a, const float* b, float* result) {
    float32x2_t va = vld1_f32(a);
    float32x2_t vb = vld1_f32(b);
    va = vsub_f32(va, vb);
    vst1_f32(result, va);
}

template<>
inline void sub<1>(const float* a, const float* b, float* result) {
    result[0] = a[0] - b[0];
}

template<int S>
inline float32x4_t l2f(const float* a, const float* b) {
    float32x4_t va = l2f<4>(a, b);
    float32x4_t vb = l2f<S - 4>(a + 4, b + 4);
    return vaddq_f32(va, vb);
}

template<>
inline float32x4_t l2f<4>(const float* a, const float* b) {
    float32x4_t va = vld1q_f32(a);
    float32x4_t vb = vld1q_f32(b);
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

template<>
inline float32x4_t l2f<3>(const float* a, const float* b) {
    float32x4_t va = {a[0], a[1], a[2], 0.0f};
    float32x4_t vb = {b[0], b[1], b[2], 0.0f};
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

template<>
inline float32x4_t l2f<2>(const float* a, const float* b) {
    float32x4_t va = {a[0], a[1], 0.0f, 0.0f};
    float32x4_t vb = {b[0], b[1], 0.0f, 0.0f};
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

template<>
inline float32x4_t l2f<1>(const float* a, const float* b) {
    float32x4_t va = {a[0], 0.0f, 0.0f, 0.0f};
    float32x4_t vb = {b[0], 0.0f, 0.0f, 0.0f};
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

template<int D>
float vector_dis(const float* a, const float* b) {
    float32x4_t x = l2f<D>(a, b);
    float32x2_t r = vadd_f32(vget_high_f32(x), vget_low_f32(x));
    return vget_lane_f32(r, 0) + vget_lane_f32(r, 1);
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

#endif //INSARM_H
