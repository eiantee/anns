//
// Created by caiqibang 00906710 on 2024/9/26.
//
#pragma once

#ifndef CONFIG_H
#define CONFIG_H

#ifdef __aarch64__
#define CONTEST 1
#endif

// #define CONTEST 1

#ifdef CONTEST
#include "insarm.h"
#define CACHE_LINE 128
#define THREADS 16
#define TEST_DIM 512
#else
#include "insx86.h"
#define CACHE_LINE 64
#define THREADS 16
#define TEST_DIM 960
#endif
#include <pthread.h>

static int thread_local my_thread_id = 0;


inline void prefetchGlobalH(const void *data, int size) {
    for (int i = 0; i < size / CACHE_LINE + 1; i++) {
        prefetch3h((char*)data + i * CACHE_LINE);
    }
}

inline void prefetchGlobalL(const void *data, int size) {
    for (int i = 0; i < size / CACHE_LINE + 1; i++) {
        prefetch3l((char*)data + i * CACHE_LINE);
    }
}

inline void prefetchLocalH(const void *data, int size) {
    for (int i = 0; i < size / CACHE_LINE + 1; i++) {
        prefetch1h((char*)data + i * CACHE_LINE);
    }
}

inline void prefetchLocalL(const void *data, int size) {
    for (int i = 0; i < size / CACHE_LINE + 1; i++) {
        prefetch1l((char*)data + i * CACHE_LINE);
    }
}

static inline void setAffinity(int thread_id) {
    pthread_t thread = pthread_self();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(thread_id, &cpuset);
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
}

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



#endif //CONFIG_H
