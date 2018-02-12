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

def main():
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

main()
