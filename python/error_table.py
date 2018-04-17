##################################
# Copyright (C) 2018 Otmar Ertl. #
# All rights reserved.           #
##################################

from math import sqrt, isnan
import csv

dataFile = '../data/error_test.csv'

algorithmDescriptionIdx = 1

def readData():
    result = []
    headers = []
    with open(dataFile , 'r') as file:
        reader = csv.reader(file, skipinitialspace=True, delimiter=';')
        rowCounter = 0
        for r in reader:
            if rowCounter <= 1:
                headers += r
            elif rowCounter % 2 == 0:
                x = []
                x += r
            else:
                x.append([int(y) for y in r])
                result.append(x)
            rowCounter += 1

    return headers,preprocess(result)

def preprocess(data):
    resultOther = []
    resultBagMinHash1Float = []
    resultBagMinHash2Float = []
    resultBagMinHash1Binary = []
    resultBagMinHash2Binary = []
    for d2 in data:
        d = d2
        if d[algorithmDescriptionIdx] == "BagMinHash1 (float)":
            d[algorithmDescriptionIdx] = "BagMinHash (float)"
            resultBagMinHash1Float.append(d)
        elif d[algorithmDescriptionIdx] == "BagMinHash2 (float)":
            d[algorithmDescriptionIdx] = "BagMinHash (float)"
            resultBagMinHash2Float.append(d)
        elif d[algorithmDescriptionIdx] == "BagMinHash1 (binary)":
            d[algorithmDescriptionIdx] = "BagMinHash (binary)"
            resultBagMinHash1Binary.append(d)
        elif d[algorithmDescriptionIdx] == "BagMinHash2 (binary)":
            d[algorithmDescriptionIdx] = "BagMinHash (binary)"
            resultBagMinHash2Binary.append(d)
        else:
            resultOther.append(d2)

    assert(len(resultBagMinHash1Float) == len(resultBagMinHash2Float))
    for i in range(0, len(resultBagMinHash1Float)):
        assert(len(resultBagMinHash1Float[i]) == len(resultBagMinHash2Float[i]))
        for j in range(0, len(resultBagMinHash1Float[i])):
            assert(len(resultBagMinHash1Float[i][j]) == len(resultBagMinHash2Float[i][j]))

    assert(len(resultBagMinHash1Binary) == len(resultBagMinHash2Binary))
    for i in range(0, len(resultBagMinHash1Binary)):
        assert(len(resultBagMinHash1Binary[i]) == len(resultBagMinHash2Binary[i]))
        for j in range(0, len(resultBagMinHash1Binary[i])):
            assert(len(resultBagMinHash1Binary[i][j]) == len(resultBagMinHash2Binary[i][j]))

    return resultOther + resultBagMinHash1Float + resultBagMinHash1Binary



headers, data = readData()

caseDescriptionIdx = 0
algorithmDescriptionIdx = 1
numIterationsIdx = 2
hashSizeIdx = 3
trueJaccardIndexIdx = 4
zeroBitHashingSimilarityIdx = 5
histogramDataIdx = 6

assert(headers[caseDescriptionIdx] == "caseDescription")
assert(headers[algorithmDescriptionIdx] == "algorithmDescription")
assert(headers[numIterationsIdx] == "numIterations")
assert(headers[hashSizeIdx] == "hashSize")
assert(headers[trueJaccardIndexIdx] == "trueJaccardIndex")
assert(headers[histogramDataIdx] == "histogramEqualSignatureComponents")
assert(headers[zeroBitHashingSimilarityIdx] == "zeroBitHashingSimilarity")

def extractCaseDescriptions(data):
    result = []
    for d in data:
        item = d[caseDescriptionIdx]
        if item not in result:
            result.append(item)
    return result

def getTrueJaccardIndex(caseDescription, data):
    for d in data:
        if d[caseDescriptionIdx] == caseDescription:
            return float(d[trueJaccardIndexIdx])

def getHistogram(caseDescription, algorithmDescription, data):
    for d in data:
        if d[caseDescriptionIdx] == caseDescription and int(d[hashSizeIdx]) == m and d[algorithmDescriptionIdx] == algorithmDescription:
            return d[histogramDataIdx]

def getEmpiricalMSE(caseDescription, m, algorithmDescription, data):
    histo = getHistogram(caseDescription, algorithmDescription, data)
    if histo is None:
        return float('nan')
    assert(m + 1 == len(histo))
    J = getTrueJaccardIndex(caseDescription, data)
    s = 0
    for k in range(0, m + 1):
        s += histo[k] * pow(k / m  - J, 2)
    return s/getN(data)

def getN(data):
    n = None
    for d in data:
        if n is None:
            n = int(d[numIterationsIdx])
        else:
            assert(n == int(d[numIterationsIdx]))
    return n

def calculateZScore(empiricalMSE, J, c, m):
    expectedMSE = J * (1 - J) / m
    expectedVarianceEmpiricalMSE = pow(expectedMSE, 2) / c * (2. - 6. / m) + expectedMSE / (c * pow(m, 2.))
    zScoreMSE = (empiricalMSE - expectedMSE) / sqrt(expectedVarianceEmpiricalMSE)
    return zScoreMSE

case_descriptions = extractCaseDescriptions(data)

m_values = [4, 16, 64, 256, 1024]

algorithms = [
    "BagMinHash (float)",
    "BagMinHash (binary)",
    "ICWS",
    "0-Bit",
    "CCWS",
    "PCWS",
    "I2CWS"
]

algorithm_labels = {
    "BagMinHash (float)" : "BagMinHash (float)",
    "BagMinHash (binary)" : "BagMinHash (binary)",
    "ICWS" : "\\acs*{ICWS} \\cite{Ioffe2010}",
    "I2CWS" : "\\acs*{I2CWS} \\cite{Wu2017}",
    "0-Bit" : "0-bit \\cite{Li2015}",
    "PCWS" : "\\acs*{PCWS} \\cite{Wu2017a}",
    "CCWS" : "\\acs*{CCWS} \\cite{Wu2016}"
}

redLimit = 3.

print("\\begin{tabular}{lrr" + (2*len(algorithms))*"r" + "}")
print("\\toprule")
print("& &")
for alg in algorithms:
    print("& \\multicolumn{2}{c}{" + algorithm_labels[alg] + "}")
print("\\\\")
i = 4
for alg in algorithms:
    print("\\cmidrule(l){" + str(i) + "-" + str(i+1) + "}")
    i += 2

print("test case & \\symHashSize & $\\symExpectation(\\symEmpiricalMSE)$")
for alg in algorithms:
    print("& $\\symEmpiricalMSE$ & $\\symZScore$-score")
print("\\\\")

n = getN(data)

for case_description in case_descriptions:
    print("\\midrule")
    i = 0
    for m in m_values:
        J = getTrueJaccardIndex(case_description, data)
        if i == 0:
            print("\\multirowcell{4}[1em][l]{" + case_description + " \\\\ " + "$\\symJaccard = " + "\\num[group-digits = false]{" + "{:.6g}".format(J) + "}" + "$}")
        i += 1
        print("& " + str(m))

        expectedMSE = J*(1.-J)/m
        print("& \\numsci{" + ' {:.2E}'.format(expectedMSE) + "}")

        for alg in algorithms:

            mse = getEmpiricalMSE(case_description, m, alg, data)

            z = calculateZScore(mse, J, n, m)

            print("&")
            if not isnan(mse):
                print("\\numsci{" + ' {:.2E}'.format(mse) + "}")
            else:
                print("N/A")

            print("&")
            if not isnan(z):
                if (abs(z) >= redLimit):
                    print("\\leavevmode\\color{red}\\bf")
                    if (abs(z) >= 10):
                        print("\\numsci{" + ' {:.2E}'.format(z) + "}")
                    else:
                        print("\\num{" + ' {:.2f}'.format(z) + "}")
                else:
                    print("\\num{" + ' {:.2f}'.format(z) + "}")
            else:
                print("N/A")


        print("\\\\")

print("\\bottomrule")
print("\\end{tabular}")
