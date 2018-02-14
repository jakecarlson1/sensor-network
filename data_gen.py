import time
import matplotlib.pyplot as plt
from Processing.objects.topology import *

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
    # times_brute, avg_degs_brute = varryAvgDeg(method="brute")

    plt.plot(avg_degs_cell, times_cell, 'r-', label="Cell")
    plt.plot(avg_degs_sweep, times_sweep, 'b-', label="Sweep")
    # plt.plot(avg_degs_brute, times_brute, 'g-', label="Brute")
    plt.xlabel("Average Degree")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for |V| = 32,000")
    plt.legend(loc=2)
    plt.show()

def main():
    # plotForVarNodes()
    plotForVarAvgDeg()

main()
