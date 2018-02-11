##################################
# Copyright (C) 2018 Otmar Ertl. #
# All rights reserved.           #
##################################

import matplotlib.pyplot as plt
import csv
import os

dataDir = '../data/'
resultFilePrefix = 'performance_test_result_'

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='cm10')

def readData():
    result = []
    for file in os.listdir(dataDir):
        if file.startswith(resultFilePrefix):
            filename = os.path.join(dataDir, file)
            with open(filename , 'r') as file:
                reader = csv.reader(file, skipinitialspace=True, delimiter=';')
                for r in reader:
                    result.append(r)

    return result

def drawChart(data, hashSize, distribution, title):
    bag_min_hash_1_unsorted = []
    bag_min_hash_1_sorted_ascending = []
    bag_min_hash_1_sorted_descending = []
    bag_min_hash_2_unsorted = []
    bag_min_hash_2_sorted_ascending = []
    bag_min_hash_2_sorted_descending = []

    for d in data:
        if d[0] != distribution:
            continue
        if int(d[3]) != hashSize:
            continue

        dataSize = int(d[4])
        avgMaxHeapSize = float(d[6])

        if d[1] == "bag_min_hash_1":
            bag_min_hash_1_unsorted.append((dataSize, avgMaxHeapSize))
        elif d[1] == "bag_min_hash_2":
            bag_min_hash_2_unsorted.append((dataSize, avgMaxHeapSize))
        elif d[1] == "icws":
            pass
        else:
            assert(false)

    bag_min_hash_1_unsorted = sorted(bag_min_hash_1_unsorted, key=lambda l: l[0])
    bag_min_hash_2_unsorted = sorted(bag_min_hash_2_unsorted, key=lambda l: l[0])

    m = bag_min_hash_1_unsorted[-1][0]
    assert m == bag_min_hash_2_unsorted[-1][0]

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(3.75, 2.5)

    ax.set_xscale("log", basex=10)
    ax.set_yscale("log", basey=10)
    ax.set_xlim([1,1e7])
    ax.set_ylim([1e3,1e7])

    ax.set_ylabel(r"space")
    ax.set_xlabel(r"$n$")

    ax.set_title(r"$m=" + str(hashSize) + "$", fontsize=12)

    ax.plot([x[0] for x in bag_min_hash_1_unsorted], [x[1] for x in bag_min_hash_1_unsorted], marker='.', label = r"BagMinHash 1", linestyle="solid", color="#7fcdbb",markeredgecolor = 'none', zorder=10)
    ax.plot([x[0] for x in bag_min_hash_2_unsorted], [x[1] for x in bag_min_hash_2_unsorted], marker='.', label = r"BagMinHash 2", linestyle="solid", color="#2c7fb8",markeredgecolor = 'none', zorder=10)
    leg = ax.legend(loc=2,ncol=1,fontsize=10,numpoints = 1)
    leg.get_frame().set_linewidth(0.0)
    leg.get_frame().set_facecolor('none')

    fig.savefig('../paper/space_chart_' + str(hashSize) + '_' + distribution + '.svg', format='svg', dpi=1200,bbox_inches='tight',pad_inches=0)
    plt.close(fig)

distributions = {"exponential_lambda_1":r"$\text{Exponential}(\text{rate}=1)$"}
hashSizes = [256, 1024, 4096, 16384]

data = readData()

for hashSize in hashSizes:
    for distribution in distributions:
        drawChart(data, hashSize, distribution, distributions[distribution])
