//##################################
//# Copyright (C) 2018 Otmar Ertl. #
//# All rights reserved.           #
//##################################

#include "weighted_minwise_hashing.hpp"

#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;
using namespace bmh;

typedef XXHash64 hash_algorithm;

template<typename H>
void testCase(uint64_t dataSize, uint32_t hashSize, uint64_t numCycles, const vector<vector<pair<uint64_t, double>>>& testData, const H& h, const string& distributionLabel, const string& algorithmLabel) {

    chrono::steady_clock::time_point tStart = chrono::steady_clock::now();
    uint64_t sumMaxSpace = 0;
    for (const auto& data :  testData) {
        assert(data.size() == dataSize);
        WeightedHashResult result = h(data, hashSize);
        sumMaxSpace += result.maxSpace;
    }
    chrono::steady_clock::time_point tEnd = chrono::steady_clock::now();

    assert(numCycles = testData.size());

    double avgHashTime = chrono::duration_cast<chrono::duration<double>>(tEnd - tStart).count()/numCycles;
    double avgMaxSpace = sumMaxSpace / static_cast<double>(numCycles);

    cout << setprecision(numeric_limits< double >::max_digits10) << scientific;
    cout << distributionLabel << ";";
    cout << algorithmLabel << ";";
    cout << numCycles << ";";
    cout << hashSize << ";";
    cout << dataSize << ";";
    cout << avgHashTime << ";";
    cout << avgMaxSpace << endl << flush;
}

template <typename DIST, typename GEN> void testDistribution(DIST& dist, GEN& rng, const string& distributionLabel, uint32_t hashSize, uint64_t dataSize, uint64_t numCycles, bool calculateIcws) {

    // generate test data
    vector<vector<pair<uint64_t, double>>> testData(numCycles);
    for (uint64_t i = 0; i < numCycles; ++i) {

        vector<pair<uint64_t, double>> d(dataSize);
        for (uint64_t j = 0; j < dataSize; ++j) {
            uint64_t data = rng();
            double weight = dist(rng);
            assert(weight <= numeric_limits<float>::max());
            d[j] = make_pair(data, weight);
        }
        testData[i] = d;
    }

    testCase(dataSize, hashSize, numCycles, testData, bag_min_hash_1<FloatWeightDiscretization, hash_algorithm>, distributionLabel, "bag_min_hash_1");
    testCase(dataSize, hashSize, numCycles, testData, bag_min_hash_2<FloatWeightDiscretization, hash_algorithm>, distributionLabel, "bag_min_hash_2");

    if (calculateIcws) {
        testCase(dataSize, hashSize, numCycles, testData, improved_consistent_weighted_hashing<hash_algorithm>, distributionLabel, "icws");
    }
}

int main(int argc, char* argv[]) {

    uint64_t numCycles = 100;

    assert(argc==6);
    uint64_t seed = atol(argv[1]);
    uint32_t hashSize = atoi(argv[2]);
    uint64_t dataSize = atol(argv[3]);
    string distribution(argv[4]);
    bool calculateIcws(atoi(argv[5]));

    mt19937_64 rng(seed);

    if (distribution == "exponential_lambda_1") {
        exponential_distribution<double> dist;
        testDistribution(dist, rng, distribution, hashSize, dataSize, numCycles, calculateIcws);
    }
    else if (distribution == "uniform_min_0_max_1") {
        uniform_real_distribution<double> dist;
        testDistribution(dist, rng, distribution, hashSize, dataSize, numCycles, calculateIcws);
    }
    else if (distribution == "weibull_scale_1_shape_0_1") {
        weibull_distribution<double> dist(0.1,1);
        testDistribution(dist, rng, distribution, hashSize, dataSize, numCycles, calculateIcws);
    }
    else {
        assert(false);
    }

    return 0;
}
