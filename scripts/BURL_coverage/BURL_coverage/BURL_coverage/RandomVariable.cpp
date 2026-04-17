#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <thread>
#include "RandomVariable.hpp"



RandomVariable::RandomVariable(void) {

    static std::atomic<uint64_t> seed_counter{0};
    
    uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    seed += seed_counter.fetch_add(1, std::memory_order_relaxed);
    seed ^= std::hash<std::thread::id>{}(std::this_thread::get_id());
    uint32_t seed32 = static_cast<uint32_t>(seed ^ (seed >> 32));
    initialize(seed32);
}

RandomVariable::RandomVariable(uint32_t seed) {

    initialize(seed);
}

uint32_t RandomVariable::extractU32(void) {

    int i = index;
    if (index >= N)
        {
        twist();
        i = index;
        }

    uint32_t y = mt[i];
    index = i + 1;

    y ^= (mt[i] >> U);
    y ^= (y << S) & B;
    y ^= (y << T) & C;
    y ^= (y >> L);

    return y;
}


void RandomVariable::initialize(uint32_t seed) {

    mt[0] = seed;
    for (uint32_t i=1; i<N; i++)
        {
        mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
        }
    index = N;
    
    for (size_t i=0; i<10000; i++)
        extractU32();
}



void RandomVariable::twist(void) {

    for (uint32_t i=0; i<N; i++)
        {
        uint32_t x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);
        uint32_t xA = x >> 1;

        if ( x & 0x1 )
            xA ^= A;

        mt[i] = mt[(i + M) % N] ^ xA;
        }
    index = 0;
}

double RandomVariable::uniformRv(void) {

    return (double)extractU32() / UINT32_MAX;
}
