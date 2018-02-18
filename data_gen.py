import time
import matplotlib.pyplot as plt
import numpy as np
from Processing.objects.topology import *
from collections import Counter

# benchmarks (num_nodes, avg_deg)
SQUARE_BENCHMARKS = [(1000,32), (8000,64), (16000,32), (64000,64), (64000,128),
                     (128000,64),, (128000, 128)]
DISK_BENCHMARKS = [(8000,64), (64000,64), (64000,128)]
SPHERE_BENCHMARKS = [(16000,64), (32000,128), (64000,128)]

def varryNumNodes(method="cell"):
    print "-------- {} --------".format(method)
    times = []
    nodes = range(4000,68000,4000)
    for n in nodes:
        topology = Square()
        topology.num_nodes = n
        topology.avg_deg = 16

        time_0 = time.clock()

        topology.generateNodes()
        topology.findEdges(method=method)

        time_0 = time.clock() - time_0
        times.append(time_0)

        print "{}: {} s".format(n, time_0)

        if time_0 > 90:
            break

    return times, nodes[:len(times)]

def varryAvgDeg(method="cell"):
    print "-------- {} --------".format(method)
    times = []
    avg_degs = range(2,34,2)
    for d in avg_degs:
        topology = Square()
        topology.num_nodes = 32000
        topology.avg_deg = d

        time_0 = time.clock()

        topology.generateNodes()
        topology.findEdges(method=method)

        time_0 = time.clock() - time_0
        times.append(time_0)

        print "{}: {} s".format(d, time_0)

        if time_0 > 90:
            break

    return times, avg_degs[:len(times)]

def plotForVarNodes():
    times_cell, nodes_cell = varryNumNodes(method="cell")
    times_sweep, nodes_sweep = varryNumNodes(method="sweep")
    times_brute, nodes_brute = varryNumNodes(method="brute")

    plt.plot(nodes_cell, times_cell, 'r-', label="Cell")
    plt.plot(nodes_sweep, times_sweep, 'b-', label="Sweep")
    plt.plot(nodes_brute, times_brute, 'g-', label="Brute")
    plt.xlabel("Num nodes")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for Average Degree = 16")
    plt.legend(loc=2)
    plt.show()

def plotForVarAvgDeg():
    times_cell, avg_degs_cell = varryAvgDeg(method="cell")
    times_sweep, avg_degs_sweep = varryAvgDeg(method="sweep")

    plt.plot(avg_degs_cell, times_cell, 'r-', label="Cell")
    plt.plot(avg_degs_sweep, times_sweep, 'b-', label="Sweep")
    plt.xlabel("Average Degree")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for |V| = 32,000")
    plt.legend(loc=2)
    plt.show()

def plotDistributionOfDegrees():
    topology = Square()
    topology.num_nodes = 32000
    topology.avg_deg = 16

    topology.generateNodes()
    topology.findEdges(method="cell")

    c = Counter([len(v) for v in topology.edges.values()])
    labels = [x[0] for x in c.items()]
    values = [x[1] for x in c.items()]
    indexes = np.arange(len(labels))
    plt.bar(indexes, values)
    plt.xticks(indexes + 0.5, labels)
    plt.xlabel("Degree")
    plt.ylabel("Number of Occurances")
    plt.title("Distribution Of Degrees of Nodes")
    plt.show()


def main():
    # plotForVarNodes()
    # plotForVarAvgDeg()
    # plotDistributionOfDegrees()

main()
