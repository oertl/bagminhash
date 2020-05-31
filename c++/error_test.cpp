//##################################
//# Copyright (C) 2018 Otmar Ertl. #
//# All rights reserved.           #
//##################################

#include "bitstream_random.hpp"
#include "weighted_minwise_hashing.hpp"
#include "data_generation.hpp"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace bmh;

typedef XXHash64 hash_algorithm;

template<typename H>
void testCase(const string& caseDescription, const Weights& w, const H& h, const string& algorithmDescription, uint64_t m, uint64_t n) {

    auto dataSizes = w.getSize();

    double J = w.getJaccardIndex();

    uint64_t seedSize = 256;

    // values from random.org
    seed_seq initialSeedSequence{
        UINT32_C(0xf12873eb), UINT32_C(0x14b77d0b), UINT32_C(0x5bf770b4), UINT32_C(0x51d8975e),
        UINT32_C(0xe6dc72f6), UINT32_C(0x5f64d296), UINT32_C(0xe2946980), UINT32_C(0x8eb55f99)};

    mt19937 initialRng(initialSeedSequence);
    vector<uint32_t> seeds(n*seedSize);
    generate(seeds.begin(), seeds.end(), initialRng);

    vector<uint32_t> numEquals(m+1);

    #pragma omp parallel for
    for (uint64_t i = 0; i < n; ++i) {

        seed_seq seedSequence(seeds.begin() + i * seedSize, seeds.begin() + (i + 1) * seedSize);
        mt19937_64 rng(seedSequence);

        const tuple<vector<tuple<uint64_t, double>>,vector<tuple<uint64_t, double>>> data = generateData(rng, w);

        const vector<tuple<uint64_t, double>>& d1 = get<0>(data);
        const vector<tuple<uint64_t, double>>& d2 = get<1>(data);

        assert(get<0>(dataSizes) == d1.size());
        assert(get<1>(dataSizes) == d2.size());

        WeightedHashResult h1 = h(d1, m);
        WeightedHashResult h2 = h(d2, m);

        uint32_t numEqual = 0;
        for (uint32_t j = 0; j < m; ++j)  {
            if (h1.hashValues[j] == h2.hashValues[j]) {
                numEqual += 1;
            }
        }

        #pragma omp atomic
        numEquals[numEqual] += 1;
    }

    cout << setprecision(numeric_limits< double >::max_digits10) << scientific;
    cout << caseDescription << ";";
    cout << algorithmDescription << ";";
    cout << n << ";";
    cout << m << ";";
    cout << J << endl;
    bool first = true;
    for(uint32_t x : numEquals) {
        if (!first) {
            cout << ";";
        }
        else {
            first = false;
        }
        cout << x;
    }
    cout << endl;
    cout << flush;
}

void testCase(const string& caseDescription, const Weights& w, uint32_t hashSize, uint64_t numIterations) {
    testCase(caseDescription, w, bag_min_hash_1<FloatWeightDiscretization, hash_algorithm>, "BagMinHash1 (float)", hashSize, numIterations);
    testCase(caseDescription, w, bag_min_hash_2<FloatWeightDiscretization, hash_algorithm>, "BagMinHash2 (float)", hashSize, numIterations);
    if (w.allWeightsZeroOrOne()) {
        testCase(caseDescription, w, bag_min_hash_1<BinaryWeightDiscretization, hash_algorithm>, "BagMinHash1 (binary)", hashSize, numIterations);
        testCase(caseDescription, w, bag_min_hash_2<BinaryWeightDiscretization, hash_algorithm>, "BagMinHash2 (binary)", hashSize, numIterations);
    }
    testCase(caseDescription, w, improved_consistent_weighted_hashing<hash_algorithm>, "ICWS", hashSize, numIterations);
    testCase(caseDescription, w, zero_bit_consistent_weighted_sampling<hash_algorithm>, "0-Bit", hashSize, numIterations);
    testCase(caseDescription, w, improved_squared_consistent_weighted_hashing<hash_algorithm>, "I2CWS", hashSize, numIterations);
    testCase(caseDescription, w, practical_consistent_weighted_hashing<hash_algorithm>, "PCWS", hashSize, numIterations);
    testCase(caseDescription, w, canonical_consistent_weighted_hashing<hash_algorithm>, "CCWS", hashSize, numIterations);
}

int main(int argc, char* argv[]) {
    cout << "caseDescription" << ";";
    cout << "algorithmDescription" << ";";
    cout << "numIterations" << ";";
    cout << "hashSize" << ";";
    cout << "trueJaccardIndex";
    cout << endl;
    cout << "histogramEqualSignatureComponents";
    cout << endl;
    cout << flush;

    uint32_t hashSizes[] = {4, 16, 64, 256, 1024, 4096};

    uint64_t numIterations = 10000;

    vector<Weights> cases = {getWeightsCase1(), getWeightsCase2(), getWeightsCase3(), getWeightsCase4(), getWeightsCase5(), getWeightsCase6(), getWeightsCase7(), getWeightsCase8(), getWeightsCase9()};

    for (uint32_t hashSize : hashSizes) {
        for (const Weights w : cases) {
            testCase(w.getLatexDescription(), w, hashSize, numIterations);
        }
    }

    return 0;
}
