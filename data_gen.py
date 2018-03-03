import time
import matplotlib.pyplot as plt
import numpy as np
import pickle
from Processing.objects.topology import *
from collections import Counter

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

def runForVarNodes():
    times_cell, nodes_cell = varryNumNodes(method="cell")
    with open('./report/data/times_cell.pkl', 'w') as f:
        pickle.dump(times_cell, f)
    with open('./report/data/nodes_cell.pkl', 'w') as f:
        pickle.dump(nodes_cell, f)

    times_sweep, nodes_sweep = varryNumNodes(method="sweep")
    with open('./report/data/times_sweep.pkl', 'w') as f:
        pickle.dump(times_sweep, f)
    with open('./report/data/nodes_sweep.pkl', 'w') as f:
        pickle.dump(nodes_sweep, f)

    times_brute, nodes_brute = varryNumNodes(method="brute")
    with open('./report/data/times_brute.pkl', 'w') as f:
        pickle.dump(times_brute, f)
    with open('./report/data/nodes_brute.pkl', 'w') as f:
        pickle.dump(nodes_brute, f)

def plotForVarNodes():
    with open('./report/data/times_cell.pkl', 'r') as f:
        times_cell = pickle.load(f)
    with open('./report/data/nodes_cell.pkl', 'r') as f:
        nodes_cell = pickle.load(f)

    with open('./report/data/times_sweep.pkl', 'r') as f:
        times_sweep = pickle.load(f)
    with open('./report/data/nodes_sweep.pkl', 'r') as f:
        nodes_sweep = pickle.load(f)

    with open('./report/data/times_brute.pkl', 'r') as f:
        times_brute = pickle.load(f)
    with open('./report/data/nodes_brute.pkl', 'r') as f:
        nodes_brute = pickle.load(f)

    plt.plot(nodes_cell, times_cell, 'r-', label="Cell")
    plt.plot(nodes_sweep, times_sweep, 'b-', label="Sweep")
    plt.plot(nodes_brute, times_brute, 'g-', label="Brute")
    plt.xlabel("Num nodes")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for Average Degree = 16")
    plt.legend(loc=2)
    plt.show()

def runForVarAvgDeg():
    times_cell, avg_degs_cell = varryAvgDeg(method="cell")
    with open('./report/data/times_cell_deg.pkl', 'w') as f:
        pickle.dump(times_cell, f)
    with open('./report/data/avg_degs_cell.pkl', 'w') as f:
        pickle.dump(avg_degs_cell, f)

    times_sweep, avg_degs_sweep = varryAvgDeg(method="sweep")
    with open('./report/data/times_sweep_deg.pkl', 'w') as f:
        pickle.dump(times_sweep, f)
    with open('./report/data/avg_degs_sweep.pkl', 'w') as f:
        pickle.dump(avg_degs_sweep, f)

def plotForVarAvgDeg():
    with open('./report/data/times_cell_deg.pkl', 'r') as f:
        times_cell = pickle.load(f)
    with open('./report/data/avg_degs_cell.pkl', 'r') as f:
        avg_degs_cell = pickle.load(f)

    with open('./report/data/times_sweep_deg.pkl', 'r') as f:
        times_sweep = pickle.load(f)
    with open('./report/data/avg_degs_sweep.pkl', 'r') as f:
        avg_degs_sweep = pickle.load(f)

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
    plt.title("Distribution Of Degrees of Nodes, Square")
    plt.show()


def main():
    # runForVarNodes()
    runForVarAvgDeg()
    # plotForVarNodes()
    # plotForVarAvgDeg()
    # plotDistributionOfDegrees()

main()
