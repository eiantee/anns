//
// Created by caiqibang 00906710 on 2024/9/26.
//
#pragma once

#ifndef CONFIG_H
#define CONFIG_H


//#define CONTEST 1

#ifdef CONTEST
#include "insarm.h"
#define CACHE_LINE 128
#else
#include "insx86.h"
#define CACHE_LINE 64
#endif
#include <pthread.h>


template<int D>
void prefetch(const void* data) {

}

static inline void setAffinity(int thread_id) {
    pthread_t thread = pthread_self();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(thread_id, &cpuset);
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
}



#endif //CONFIG_H
