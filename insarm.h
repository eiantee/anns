//
// Created by caiqibang 00906710 on 2024/9/26.
//

#ifndef INSARM_H
#define INSARM_H
#include <arm_neon.h>

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
inline void sub(float* a, float* b, float* result) {
    sub<4>(a, b, result);
    sub<S - 4>(a + 4, b + 4, result + 4);
}

inline void sub<4>(float* a, float* b, float* result) {
    float32x4_t va = vld1q_f32(a);
    float32x4_t vb = vld1q_f32(b);
    va = vsubq_f32(va, vb);
    vst1q_f32(result, va);
}

inline void sub<3>(float* a, float* b, float* result) {
    float32x4_t va = {a[0], a[1], a[2], 0.0f};
    float32x4_t vb = {b[0], b[1], b[2], 0.0f};
    va = vsubq_f32(va, vb);
    vst1q_lane_f32(result, va, 0);
    vst1q_lane_f32(result + 1, va, 1);
    vst1q_lane_f32(result + 2, va, 2);
}

inline void sub<2>(float* a, float* b, float* result) {
    float32x2_t va = vld1_f32(a);
    float32x2_t vb = vld1_f32(b);
    va = vsub_f32(va, vb);
    vst1_f32(result, va);
}

inline void sub<1>(float* a, float* b, float* result) {
    result[0] = a[0] - b[0];
}

template<int S>
inline float32x4_t l2f(float* a, float* b) {
    float32x4_t va = l2f<4>(a, b);
    float32x4_t vb = l2f<S - 4>(a + 4, b + 4);
    return vaddq_f32(va, vb);
}

inline float32x4_t l2f<4>(float* a, float* b) {
    float32x4_t va = vld1q_f32(a);
    float32x4_t vb = vld1q_f32(b);
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

inline float32x4_t l2f<3>(float* a, float* b) {
    float32x4_t va = {a[0], a[1], a[2], 0.0f};
    float32x4_t vb = {b[0], b[1], b[2], 0.0f};
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

inline float32x4_t l2f<2>(float* a, float* b) {
    float32x4_t va = {a[0], a[1], 0.0f, 0.0f};
    float32x4_t vb = {b[0], b[1], 0.0f, 0.0f};
    va = vsubq_f32(va, vb);
    return vmulq_f32(va, va);
}

inline float32x4_t l2f<1>(float* a, float* b) {
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

#endif //INSARM_H
