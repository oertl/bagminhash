//##################################
//# Copyright (C) 2018 Otmar Ertl. #
//# All rights reserved.           #
//##################################

#ifndef _DATA_GENERATION_HPP_
#define _DATA_GENERATION_HPP_

#include <vector>
#include <string>
#include <random>
#include <tuple>
#include <algorithm>

namespace bmh {

class Weights {

    const std::vector<std::tuple<double, double>> weights;
    const std::string latexDescription;

public:
    double getJaccardIndex() const {
        double d = 0;
        double n = 0;
        for(const auto& w : weights) {
            d += std::max(std::get<0>(w), std::get<1>(w));
            n += std::min(std::get<0>(w), std::get<1>(w));
        }
        return n/d;
    }

    std::tuple<size_t, size_t> getSize() const {
        size_t sizeA = 0;
        size_t sizeB = 0;
        for (const auto& w : weights) {
            if(std::get<0>(w) > 0) sizeA += 1;
            if(std::get<1>(w) > 0) sizeB += 1;
        }
        return std::make_tuple(sizeA, sizeB);
    }

    const std::vector<std::tuple<double, double>> getWeights() const {
        return weights;
    }

    const std::string& getLatexDescription() const {
        return latexDescription;
    }

    bool allWeightsZeroOrOne() const {
        for(const auto& w : weights) {
            if (
                (std::get<0>(w) != 0 && std::get<0>(w)!=1) ||
                (std::get<1>(w) != 0 && std::get<1>(w)!=1)) return false;
        }
        return true;
    }

    Weights(const std::vector<std::tuple<double, double>>& w, const std::string& d) : weights(w), latexDescription(d) {}
};

std::tuple<std::vector<std::tuple<uint64_t, double>>, std::vector<std::tuple<uint64_t, double>>> generateData(std::mt19937_64& rng, const Weights& w) {

    const auto& resultSizes = w.getSize();

    std::vector<std::tuple<uint64_t, double>> valuesA;
    std::vector<std::tuple<uint64_t, double>> valuesB;

    valuesA.reserve(std::get<0>(resultSizes));
    valuesB.reserve(std::get<1>(resultSizes));

    for(const auto& x : w.getWeights()) {
        uint64_t data = rng();
        if(std::get<0>(x) > 0) valuesA.emplace_back(data, std::get<0>(x));
        if(std::get<1>(x) > 0) valuesB.emplace_back(data, std::get<1>(x));
    }

    std::shuffle(valuesA.begin(), valuesA.end(), rng);
    std::shuffle(valuesB.begin(), valuesB.end(), rng);

    return std::make_tuple(valuesA, valuesB);
}

Weights getWeightsCase1() {
    std::vector<std::tuple<double,double>> v;
    v.emplace_back(1,10);
    return Weights(v, "$\\lbrace(1,10)\\rbrace$");
}

Weights getWeightsCase2() {
    std::vector<std::tuple<double,double>> v;
    v.emplace_back(9,10);
    return Weights(v, "$\\lbrace(9,10)\\rbrace$");
}

Weights getWeightsCase3() {
    std::vector<std::tuple<double,double>> v;
    v.emplace_back(3,20);
    v.emplace_back(30,7);
    return Weights(v, "$\\lbrace(3,20),(30,7)\\rbrace$");
}

Weights getWeightsCase4() {
    std::vector<std::tuple<double,double>> v;
    v.emplace_back(0,2);
    v.emplace_back(3,4);
    v.emplace_back(6,3);
    v.emplace_back(2,4);
    return Weights(v, "$\\lbrace(0,2),(3,4),(6,3),(2,4)\\rbrace$");
}

Weights getWeightsCase5() {
    std::vector<std::tuple<double,double>> v;
    for (int i = 0; i < 15; ++i) {
        v.emplace_back(4,2);
    }
    for (int i = 0; i < 10; ++i) {
        v.emplace_back(1,4);
    }
    for (int i = 0; i < 5; ++i) {
        v.emplace_back(12,0);
    }
    return Weights(v, "$\\lbrace(4,2)^{15},(1,4)^{10},(12,0)^{5}\\rbrace$");
}

Weights getWeightsCase6() {
    std::vector<std::tuple<double,double>> v;

    for (int i = 0; i <= 1000; ++i) {
        v.emplace_back(std::pow(1.001,i), std::pow(1.002, i));
    }
    return Weights(v, "$\\bigcup_{\\symTestCaseIndex=0}^{1000} \\lbrace({1.001}^\\symTestCaseIndex, {1.002}^{\\symTestCaseIndex})\\rbrace$");
}

Weights getWeightsCase7() {
    std::vector<std::tuple<double,double>> v;
    v.emplace_back(0,1);
    v.emplace_back(1,1);
    v.emplace_back(1,0);
    return Weights(v, "$\\lbrace(0,1),(1,0),(1,1)\\rbrace$");
}

Weights getWeightsCase8() {
    std::vector<std::tuple<double,double>> v;
    for (int i = 0; i < 30; ++i) {
        v.emplace_back(0,1);
    }
    for (int i = 0; i < 10; ++i) {
        v.emplace_back(1,0);
    }
    for (int i = 0; i < 160; ++i) {
        v.emplace_back(1,1);
    }
    return Weights(v, "$\\lbrace(0,1)^{30},(1,0)^{10},(1,1)^{160}\\rbrace$");
}

Weights getWeightsCase9() {
    std::vector<std::tuple<double,double>> v;
    for (int i = 0; i < 300; ++i) {
        v.emplace_back(0,1);
    }
    for (int i = 0; i < 500; ++i) {
        v.emplace_back(1,0);
    }
    for (int i = 0; i < 1200; ++i) {
        v.emplace_back(1,1);
    }
    return Weights(v, "$\\lbrace(0,1)^{300},(1,0)^{500},(1,1)^{1200}\\rbrace$");
}

} // namespace bmh

#endif // _DATA_GENERATION_HPP_
