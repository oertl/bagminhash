//##################################
//# Copyright (C) 2018 Otmar Ertl. #
//# All rights reserved.           #
//##################################

#ifndef _WEIGHTED_MINWISE_HASHING_HPP_
#define _WEIGHTED_MINWISE_HASHING_HPP_

#include "bitstream_random.hpp"

#include <vector>
#include <algorithm>
#include <functional>

template <typename T>
class MaxValueTracker {
    const uint32_t m;
    std::vector<T> values;

public:
    MaxValueTracker(uint32_t _m, const T& infinity) : m(_m), values((_m << 1) - 1, infinity) {}

    void update(uint32_t idx, T value) {
        assert(idx < m);
        while(value < values[idx]) {
            values[idx] = value;
            idx = m + (idx >> 1);
            if (idx >= values.size()) break;
            uint32_t leftChildIdx = (idx - m) << 1;
            uint32_t rightChildIdx = leftChildIdx + 1;
            value = std::max(values[leftChildIdx], values[rightChildIdx]);
        }
    }

    const T& max() const {
        return values.back();
    }

    const T& operator[](uint32_t idx) const {
        return values[idx];
    }
};

class BinaryWeightDiscretization {
public:
    typedef uint8_t index_type;
    typedef uint8_t weight_type;

    static const index_type maxBoundIdx = 1;

    static weight_type getBound(index_type boundIdx) {
        return boundIdx;
    }
};

template<typename I, typename W>
I calculateMaxBoundIdx() {
    static_assert(sizeof(I) == sizeof(W), "index type and weight type do not have same size");
    I i = 0;
    const W w = std::numeric_limits<W>::max();
    memcpy(&i, &w, sizeof(I));
    return i;
}

template<typename W, typename I>
class WeightDiscretization {
public:
    typedef I index_type;
    typedef W weight_type;

    static const index_type maxBoundIdx;

    static weight_type getBound(index_type boundIdx) {
        W f;
        static_assert(std::numeric_limits<W>::is_iec559, "weight_type is not iec559");
        static_assert(sizeof(weight_type) == sizeof(index_type), "weight_type and index_type do not have same size");
        memcpy(&f, &boundIdx, sizeof(index_type));
        return f;
    }
};

template<typename W, typename I>
const typename WeightDiscretization<W, I>::index_type WeightDiscretization<W, I>::maxBoundIdx = calculateMaxBoundIdx<WeightDiscretization<W, I>::index_type, WeightDiscretization<W, I>::weight_type>();

typedef WeightDiscretization<float, uint32_t> FloatWeightDiscretization;

typedef WeightDiscretization<double, uint64_t> DoubleWeightDiscretization;

struct WeightedHashResult {
    std::vector<uint64_t> hashValues;
    uint64_t maxSpace;

    WeightedHashResult(uint32_t m) : hashValues(m), maxSpace(UINT64_C(0)) {}
};

template<typename... V> class ValueProvider;

template<>
class ValueProvider<> {
public:
    static size_t size() {
        return 0;
    }

    void init(char* data) const {}
};

template<typename W, typename... V> class ValueProvider<W, V...> : ValueProvider<V...> {
    const W w;
public:

    ValueProvider(const W& _w, const V& ... _v) : ValueProvider<V...>(_v...), w(_w) {}

    static size_t size() {
        return ValueProvider<V...>::size() + sizeof(W);
    }

    void init(char* data) const {
        ValueProvider<V...>::init(data);
        memcpy(&data[ValueProvider<V...>::size()], &w, sizeof(W));
    }
};

template<typename... V>
ValueProvider<V...> collectHashData(const V&... v) {
    return ValueProvider<V...>(v...);
}

template <typename D, typename H>
class PoissonProcess {

    double point;
    double weight;
    typename D::index_type weightIdxMin;
    typename D::index_type weightIdxMax;
    typename D::weight_type boundMin;
    typename D::weight_type boundMax;
    uint32_t signatureIdx;
    BitStream<H> randomBitStream;

public:

    PoissonProcess(
        double _point,
        double _weight,
        typename D::index_type _weightIdxMin,
        typename D::index_type _weightIdxMax,
        typename D::weight_type _boundMin,
        typename D::weight_type _boundMax,
        BitStream<H>&& _randomBitStream
        )
      : point(_point),
        weight(_weight),
        weightIdxMin(_weightIdxMin),
        weightIdxMax(_weightIdxMax),
        boundMin(_boundMin),
        boundMax(_boundMax),
        signatureIdx(std::numeric_limits<uint32_t>::max()),
        randomBitStream(std::move(_randomBitStream)) {}


    PoissonProcess(BitStream<H>&& _randomBitStream, double _weight) :
        PoissonProcess(0., _weight, 0, D::maxBoundIdx, 0, D::getBound(D::maxBoundIdx), std::move(_randomBitStream)) {}

    bool splittable() const {
        return weightIdxMax > weightIdxMin + 1;
    }

    bool partiallyRelevant() const {
        return D::getBound(weightIdxMin + 1) <= weight;
    }

    bool fullyRelevant() const {
        return boundMax <= weight;
    }

    uint32_t getIndex() const {
        return signatureIdx;
    }

    double getPoint() const {
        return point;
    }

    void next(uint32_t m) {
        point += getExponential1(randomBitStream) / (static_cast<double>(boundMax) - static_cast<double>(boundMin));
        signatureIdx = getUniform(m, randomBitStream);
    }

    std::unique_ptr<PoissonProcess> split() {

        typename D::index_type weightIdxMid = (weightIdxMin + weightIdxMax) >> 1;

        double boundMid = D::getBound(weightIdxMid);

        bool inheritToLeft = getBernoulli((boundMid - static_cast<double>(boundMin)) / (static_cast<double>(boundMax) - static_cast<double>(boundMin)), randomBitStream);

        std::unique_ptr<PoissonProcess> pPrime;

        BitStream<H> bitStream(collectHashData(weightIdxMid, point), UINT64_C(0x4b06d55ba29b0826)); // constant from random.org

        if (inheritToLeft) {
            pPrime = std::make_unique<PoissonProcess>(point, weight, weightIdxMid, weightIdxMax, boundMid, boundMax, std::move(bitStream));
            weightIdxMax = weightIdxMid;
            boundMax = boundMid;
        }
        else {
            pPrime = std::make_unique<PoissonProcess>(point, weight, weightIdxMin, weightIdxMid, boundMin, boundMid, std::move(bitStream));
            weightIdxMin = weightIdxMid;
            boundMin = boundMid;
        }
        return pPrime;
    }
};

template<typename D, typename H>
struct CmpPoissonProcessPtrs
{
    bool operator()(const std::unique_ptr<PoissonProcess<D,H>>& lhs, const std::unique_ptr<PoissonProcess<D,H>>& rhs) const
    {
        return rhs->getPoint() < lhs->getPoint();
    }
};

template<typename D, typename H>
struct CmpPoissonProcessPtrsInverse
{
    bool operator()(const std::unique_ptr<PoissonProcess<D,H>>& lhs, const std::unique_ptr<PoissonProcess<D,H>>& rhs) const
    {
        return rhs->getPoint() > lhs->getPoint();
    }
};

template<typename D, typename H>
void pushHeap(std::unique_ptr<PoissonProcess<D, H>>& p, std::vector<std::unique_ptr<PoissonProcess<D,H>>>& heap, uint64_t& maxHeapSize, size_t offset = 0) {
    heap.emplace_back(std::move(p));
    std::push_heap(heap.begin() + offset, heap.end(), CmpPoissonProcessPtrs<D,H>());
    if (heap.size() > maxHeapSize) maxHeapSize = heap.size();
}

template<typename D, typename H>
std::unique_ptr<PoissonProcess<D, H>> popHeap(std::vector<std::unique_ptr<PoissonProcess<D,H>>>& heap, size_t offset = 0) {
    std::pop_heap(heap.begin() + offset, heap.end(), CmpPoissonProcessPtrs<D,H>());
    std::unique_ptr<PoissonProcess<D,H>> p = std::move(heap.back());
    heap.pop_back();
    return p;
}

static const uint64_t bagMinHashSeedA = UINT64_C(0xf331e07615a87fd7); // constant from random.org
static const uint64_t bagMinHashSeedB = UINT64_C(0xe224afad0d89c684); // constant from random.org

template<typename D, typename H>
WeightedHashResult bag_min_hash_1(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {
    assert(D::getBound(0) == 0);

    const uint8_t b = 64; // constant for b-bit minwise hashing

    std::vector<std::unique_ptr<PoissonProcess<D, H>>> heap;

    MaxValueTracker<double> h(m, std::numeric_limits<double>::infinity());
    WeightedHashResult result(m);

    for(const auto& item : data) {

        const double w = std::get<1>(item);
        if (w < D::getBound(1)) continue;

        const uint64_t d = std::get<0>(item);

        BitStream<H> bitStream(collectHashData(d), bagMinHashSeedA);
        std::unique_ptr<PoissonProcess<D,H>> p = std::make_unique<PoissonProcess<D,H>>(std::move(bitStream), w);

        p->next(m);
        if (p->fullyRelevant()) h.update(p->getIndex(), p->getPoint());

        while(p->getPoint() <= h.max()) {
            while(p->splittable() && p->partiallyRelevant()) {

                std::unique_ptr<PoissonProcess<D,H>> pPrime = p->split();

                if (p->fullyRelevant()) h.update(p->getIndex(), p->getPoint());

                if (pPrime->partiallyRelevant()) {
                    pPrime->next(m);
                    if (pPrime->fullyRelevant()) h.update(pPrime->getIndex(), pPrime->getPoint());
                    if (pPrime->getPoint() <= h.max()) pushHeap(pPrime, heap, result.maxSpace);
                }
            }

            if (p->fullyRelevant()) {
                p->next(m);
                h.update(p->getIndex(), p->getPoint());
                if (p->getPoint() <= h.max()) pushHeap(p, heap, result.maxSpace);
            }
            if (heap.empty()) break;
            p = popHeap(heap);
        }

        heap.clear();
    }

    for (uint32_t k = 0; k < m; ++k) {
        BitStream<H> bitstream(collectHashData(h[k]), bagMinHashSeedB);
        result.hashValues[k] = getUniformPow2(b, bitstream);
    }

    return result;
}

template<typename D, typename H>
WeightedHashResult bag_min_hash_2(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {
    assert(D::getBound(0) == 0);

    const uint8_t b = 64; // constant for b-bit minwise hashing

    std::vector<std::unique_ptr<PoissonProcess<D,H>>> temp;
    WeightedHashResult result(m);

    MaxValueTracker<double> h(m, std::numeric_limits<double>::infinity());

    for(const auto& item : data) {

        const double w = std::get<1>(item);
        if (w < D::getBound(1)) continue;

        const uint64_t d = std::get<0>(item);

        const size_t tempHeapOffset = temp.size();

        BitStream<H> bitStream(collectHashData(d), bagMinHashSeedA);
        std::unique_ptr<PoissonProcess<D,H>> p = std::make_unique<PoissonProcess<D,H>>(std::move(bitStream), w);

        p->next(m);
        if (p->fullyRelevant()) h.update(p->getIndex(), p->getPoint());

        while(p->getPoint() <= h.max()) {
            while(p->splittable() && p->partiallyRelevant() && !p->fullyRelevant()) {

                std::unique_ptr<PoissonProcess<D,H>> pPrime = p->split();

                if (p->fullyRelevant()) h.update(p->getIndex(), p->getPoint());

                if (pPrime->partiallyRelevant()) {
                    pPrime->next(m);
                    if (pPrime->fullyRelevant()) h.update(pPrime->getIndex(), pPrime->getPoint());
                    if (pPrime->getPoint() <= h.max()) pushHeap(pPrime, temp, result.maxSpace, tempHeapOffset);
                }
            }

            if (p->fullyRelevant()) {
                assert(p->getPoint() <= h.max());
                pushHeap(p, temp, result.maxSpace, tempHeapOffset);
                break;
            }
            if (temp.size() == tempHeapOffset) break;
            p = popHeap(temp, tempHeapOffset);
        }

        auto bufferEndIt = temp.begin() + tempHeapOffset;
        while(bufferEndIt != temp.begin() && temp.front()->getPoint() > h.max()) {
            std::pop_heap(temp.begin(), bufferEndIt, CmpPoissonProcessPtrsInverse<D,H>());
            --bufferEndIt;
        }

        for(auto heapIt = temp.begin() + tempHeapOffset; heapIt != temp.end(); ++heapIt) {
            if ((*heapIt)->getPoint() <= h.max()) {
                *bufferEndIt = std::move(*heapIt);
                ++bufferEndIt;
                std::push_heap(temp.begin(), bufferEndIt, CmpPoissonProcessPtrsInverse<D,H>());
            }
        }
        temp.erase(bufferEndIt, temp.end());
    }

    std::make_heap(temp.begin(), temp.end(), CmpPoissonProcessPtrs<D,H>());

    while(!temp.empty()) {

        std::unique_ptr<PoissonProcess<D,H>> p = popHeap(temp);
        if (p->getPoint() > h.max()) break;

        while(p->splittable() && p->partiallyRelevant()) {

            std::unique_ptr<PoissonProcess<D,H>> pPrime = p->split();

            if (p->fullyRelevant()) h.update(p->getIndex(), p->getPoint());

            if (pPrime->partiallyRelevant()) {
                pPrime->next(m);
                if (pPrime->fullyRelevant()) h.update(pPrime->getIndex(), pPrime->getPoint());
                if (pPrime->getPoint() <= h.max()) pushHeap(pPrime, temp, result.maxSpace);
            }
        }

        if (p->fullyRelevant()) {
            p->next(m);
            h.update(p->getIndex(), p->getPoint());
            if (p->getPoint() <= h.max()) pushHeap(p, temp, result.maxSpace);
        }
    }

    for (uint32_t k = 0; k < m; ++k) {
        BitStream<H> bitstream(collectHashData(h[k]), bagMinHashSeedB);
        result.hashValues[k] = getUniformPow2(b, bitstream);
    }

    return result;
}

// see Ioffe, Sergey. "Improved consistent sampling, weighted minhash and l1 sketching." Data Mining (ICDM), 2010 IEEE 10th International Conference on. IEEE, 2010.
template <typename H>
WeightedHashResult improved_consistent_weighted_hashing(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {
    const uint8_t b = 64; // constant for b-bit minwise hashing

    std::vector<double> aVec(m, std::numeric_limits<double>::infinity());
    std::vector<uint64_t> dVec(m);
    std::vector<uint64_t> yVec(m);

    WeightedHashResult result(m);

    for(const auto& item : data) {

        const uint64_t d = std::get<0>(item);
        const double s = std::get<1>(item);

        if (s == 0) continue;

        const double logS = std::log(s);

        BitStream<H> bitstream(collectHashData(d), UINT64_C(0x87609608d2a48b5d)); // constant from random.org

        for (uint32_t k = 0; k < m; ++k) {

            double r = getGamma21(bitstream);
            double c = getGamma21(bitstream);
            double beta = getUniformDouble(bitstream);

            double t = std::floor(logS / r + beta);
            double y = std::exp(r * (t - beta));
            double a = c / (y * std::exp(r));

            if (a < aVec[k]) {
                aVec[k] = a;
                dVec[k] = d;
                yVec[k] = y;
            }
        }
    }

    for (uint32_t k = 0; k < m; ++k) {
        BitStream<H> bitstream(collectHashData(dVec[k], yVec[k]), UINT64_C(0xbf235dea3db9c393)); // constant from random.org
        result.hashValues[k] = getUniformPow2(b, bitstream);
    }

    return result;
}

// see Wu, Wei, et al. "Canonical Consistent Weighted Sampling for Real-Value Weighted Min-Hash." Data Mining (ICDM), 2016 IEEE 16th International Conference on. IEEE, 2016.
template <typename H>
WeightedHashResult canonical_consistent_weighted_hashing(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {
    const uint8_t b = 64; // constant for b-bit minwise hashing

    std::vector<double> aVec(m, std::numeric_limits<double>::infinity());
    std::vector<uint64_t> dVec(m);
    std::vector<uint64_t> yVec(m);

    WeightedHashResult result(m);

    for(const auto& item : data) {

        const uint64_t d = std::get<0>(item);
        const double s = std::get<1>(item);

        if (s == 0) continue;

        BitStream<H> bitstream(collectHashData(d), UINT64_C(0xc9116756125c6267)); // constant from random.org

        for (uint32_t k = 0; k < m; ++k) {

            double beta = getUniformDouble(bitstream);
            double r = getBeta21(bitstream);
            double c = getGamma21(bitstream);

            double t = std::floor(s / r + beta);
            double y = r * (t - beta);
            double a = c / y - 2 * r * c;

            if (a < aVec[k]) {
                aVec[k] = a;
                dVec[k] = d;
                yVec[k] = y;
            }
        }
    }

    for (uint32_t k = 0; k < m; ++k) {
        BitStream<H> bitstream(collectHashData(dVec[k], yVec[k]), UINT64_C(0xa5c48ff7b4004c41)); // constant from random.org
        result.hashValues[k] = getUniformPow2(b, bitstream);
    }

    return result;
}

// see Wu, Wei, et al. "Consistent Weighted Sampling Made More Practical." Proceedings of the 26th International Conference on World Wide Web. International World Wide Web Conferences Steering Committee, 2017.
template <typename H>
WeightedHashResult practical_consistent_weighted_hashing(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {
    const uint8_t b = 64; // constant for b-bit minwise hashing

    std::vector<double> aVec(m, std::numeric_limits<double>::infinity());
    std::vector<uint64_t> dVec(m);
    std::vector<uint64_t> yVec(m);

    WeightedHashResult result(m);

    for(const auto& item : data) {

        const uint64_t d = std::get<0>(item);
        const double s = std::get<1>(item);

        if (s == 0) continue;

        const double logS = std::log(s);

        BitStream<H> bitstream(collectHashData(d), UINT64_C(0xbe46368ee398beee)); // constant from random.org

        for (uint32_t k = 0; k < m; ++k) {

            double u1 = getUniformDouble(bitstream);
            double u2 = getUniformDouble(bitstream);
            double beta = getUniformDouble(bitstream);
            double x = getUniformDouble(bitstream);

            double gamma = -std::log(u1 * u2);
            double t = std::floor(logS / gamma + beta);
            double y = std::exp(gamma * (t - beta));
            double a = -std::log(x) / (y / u1);

            if (a < aVec[k]) {
                aVec[k] = a;
                dVec[k] = d;
                yVec[k] = y;
            }
        }
    }

    for (uint32_t k = 0; k < m; ++k) {
        BitStream<H> bitstream(collectHashData(dVec[k], yVec[k]), UINT64_C(0x50da48973b000da9)); // constant from random.org
        result.hashValues[k] = getUniformPow2(b, bitstream);
    }

    return result;
}

// see Wu, Wei, et al. "Improved Consistent Weighted Sampling Revisited." arXiv preprint arXiv:1706.01172 (2017).
template <typename H>
WeightedHashResult improved_squared_consistent_weighted_hashing(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {
    const uint8_t b = 64; // constant for b-bit minwise hashing

    std::vector<double> aVec(m, std::numeric_limits<double>::infinity());
    std::vector<uint64_t> dVec(m);
    std::vector<uint64_t> yVec(m);

    WeightedHashResult result(m);

    for(const auto& item : data) {

        const uint64_t d = std::get<0>(item);
        const double s = std::get<1>(item);

        if (s == 0) continue;

        const double logS = std::log(s);

        BitStream<H> bitstream(collectHashData(d), UINT64_C(0xb30eb19e5e572b46)); // constant from random.org

        for (uint32_t k = 0; k < m; ++k) {

            double r1 = getGamma21(bitstream);
            double r2 = getGamma21(bitstream);
            double beta1 = getUniformDouble(bitstream);
            double beta2 = getUniformDouble(bitstream);
            double c = getGamma21(bitstream);

            double t2 = std::floor(logS / r2 + beta2);
            double z = std::exp(r2 * (t2 - beta2 + 1));
            double a = c / z;

            if (a < aVec[k]) {
                aVec[k] = a;
                double t1 = std::floor(logS / r1 + beta1);
                double y = std::exp(r1 * (t1 - beta1));
                dVec[k] = d;
                yVec[k] = y;
            }
        }
    }

    for (uint32_t k = 0; k < m; ++k) {
        BitStream<H> bitstream(collectHashData(dVec[k], yVec[k]), UINT64_C(0xdff981675040e7bc)); // constant from random.org
        result.hashValues[k] = getUniformPow2(b, bitstream);
    }

    return result;
}

// see Li, Ping. "0-bit consistent weighted sampling." Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. ACM, 2015.
template <typename H>
WeightedHashResult zero_bit_consistent_weighted_sampling(const std::vector<std::tuple<uint64_t,double>>& data, const uint32_t m) {

    std::vector<double> aVec(m, std::numeric_limits<double>::infinity());
    std::vector<uint64_t> dVec(m);

    WeightedHashResult result(m);

    for(const auto& item : data) {

        const uint64_t d = std::get<0>(item);
        const double s = std::get<1>(item);

        if (s == 0) continue;

        const double logS = std::log(s);

        BitStream<H> bitstream(collectHashData(d), UINT64_C(0x95f3ee483861a892)); // constant from random.org

        for (uint32_t k = 0; k < m; ++k) {

            double r = getGamma21(bitstream);
            double c = getGamma21(bitstream);
            double beta = getUniformDouble(bitstream);

            double t = std::floor(logS / r + beta);
            double y = std::exp(r * (t - beta));
            double a = c / (y * std::exp(r));

            if (a < aVec[k]) {
                aVec[k] = a;
                dVec[k] = d;
            }
        }
    }

    result.hashValues = std::move(dVec);

    return result;
}

#endif // _WEIGHTED_MINWISE_HASHING_HPP_
