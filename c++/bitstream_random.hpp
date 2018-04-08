//##################################
//# Copyright (C) 2018 Otmar Ertl. #
//# All rights reserved.           #
//##################################

#ifndef _BIT_STREAM_RANDOM_HPP_
#define _BIT_STREAM_RANDOM_HPP_

#include "exponential_distribution.hpp"

#include "xxhash/xxhash.h"

#include <cmath>
#include <cassert>
#include <cstring>
#include <memory>

static constexpr double maxReciprocal = 1. / (UINT64_C(1) << 52);

// uniform distributed double value from (0, 1]
template<typename T> double getUniformDouble(T& bitstream) {
    return (bitstream(52) + 1) * maxReciprocal;
}

template<typename T> double getExponential1(T& bitstream) {
    return ziggurat::getExponential(bitstream);
}

template<typename T> double getGamma21(T& bitstream) {
    return getExponential1(bitstream) + getExponential1(bitstream);
}

template<typename T> double getBeta21(T& bitstream) {
    return std::sqrt(getUniformDouble(bitstream));
}

template<typename T> bool getBernoulli(double successProbability, T& bitStream) {
    while(true) {
        if (successProbability == 0) return false;
        if (successProbability == 1) return true;
        bool b = successProbability > 0.5;
        if (bitStream()) return b;
        successProbability += successProbability;
        if (b) successProbability -= 1;
    }
}

// see Lumbroso, Jeremie. "Optimal discrete uniform generation from coin flips, and applications." arXiv preprint arXiv:1304.1916 (2013).
template<typename T> uint64_t getUniform(uint64_t n, T& bitstream) {
    assert(n > 0);
    uint64_t v = 1;
    uint64_t c = 0;
    while(true) {
        v <<= 1;
        c <<= 1;
        c += bitstream();
        if (v >= n) {
            if (c < n) {
                return c;
            }
            else {
                v -= n;
                c -= n;
            }
        }
    }
}

template<typename T> uint64_t getUniformPow2(uint8_t numBits, T& bitstream) {
    return bitstream(numBits);
}

class XXHash64 {
public:
    static uint64_t calculateHash(const char* data, size_t length, uint64_t seed) {
        return XXH64(data, length, seed);
    }
};

struct BitMasks {
    constexpr BitMasks() : masks() {
        masks[0] = 0;
        for (uint8_t i = 1; i <= 63; ++i) masks[i] = (UINT64_C(1) << i) - UINT64_C(1);
        masks[64] = UINT64_C(0xFFFFFFFFFFFFFFFF);
    }

    uint64_t masks[65];
};

static constexpr BitMasks BIT_MASKS;

template<typename R>
class BitStream {

    static const uint32_t FNV_OFFSET;
    static const uint32_t FNV_PRIME;

    size_t dataSize;
    std::unique_ptr<char[]> data;
    uint64_t seed;
    uint64_t hashBits;
    uint8_t availableBits;

    void nextHash() {
        uint32_t tmp;
        memcpy(&tmp, data.get(), sizeof(uint32_t));
        tmp *= FNV_PRIME;
        memcpy(data.get(), &tmp, sizeof(uint32_t));
        hashBits = R::calculateHash(data.get(), dataSize, seed);
    }
public:

    BitStream(const BitStream& p) = delete;
    BitStream& operator=(const BitStream&) = delete;
    BitStream(BitStream&& p) = default;
    BitStream& operator=(BitStream&&) = default;

    template<typename I>
    BitStream(const I& valueProvider, uint64_t _seed) : dataSize(valueProvider.size() +  sizeof(uint32_t)), data(new char[dataSize]), seed(_seed), hashBits(0), availableBits(0) {
        memcpy(data.get(), &FNV_OFFSET, sizeof(uint32_t));
        valueProvider.init(&data[sizeof(uint32_t)]);
    }

    bool operator()() {
        if (availableBits == 0) {
            nextHash();
        }
        bool result = (hashBits & UINT64_C(1));
        hashBits >>= 1;
        availableBits -= 1;
        availableBits &= UINT8_C(0x3F);
        return result;
    }

    uint64_t operator()(uint8_t numBits) {
        assert(numBits <= 64);
        uint64_t result = 0;
        uint8_t requiredBits = numBits;
        if(numBits > availableBits) {
            result = (hashBits & BIT_MASKS.masks[availableBits]);
            result <<= (numBits - availableBits);
            nextHash();
            requiredBits -= availableBits;
        }
        result |= (hashBits & BIT_MASKS.masks[requiredBits]);
        hashBits >>= requiredBits;
        availableBits -= numBits;
        availableBits &= UINT8_C(0x3F);
        return result;
    }
};

// see https://en.wikipedia.org/wiki/Fowler-Noll-Vo_hash_function
template<typename R> const uint32_t BitStream<R>:: FNV_OFFSET = 0x811c9dc5;
template<typename R> const uint32_t BitStream<R>:: FNV_PRIME = (1 << 24) + (1 << 8) + 0x93;

#endif // _BIT_STREAM_RANDOM_HPP_
