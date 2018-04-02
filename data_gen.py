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

    c = np.polyfit(nodes_cell, times_cell, 1)
    c_n = np.poly1d(c)
    s = np.polyfit(nodes_sweep, times_sweep, 2)
    s_n = np.poly1d(s)
    b = np.polyfit(nodes_brute, times_brute, 2)
    b_n = np.poly1d(b)

    print "cell: t = %.6fn + (%.6f)"%(c[0],c[1])
    print "sweep: t = %.6fn^2 + %.6fn + (%.6f)"%(s[0],s[1],s[2])
    print "brute: t = %.6fn^2 + %.6fn + (%.6f)"%(b[0],b[1],b[2])

    plt.plot(nodes_cell, times_cell, 'r-', label="Cell")
    plt.plot(nodes_cell, c_n(nodes_cell), 'r--')
    plt.text(42000, 6, '$t_{cell} = %.6fn + (%.6f)$'%(c[0],c[1]))
    plt.plot(nodes_sweep, times_sweep, 'b-', label="Sweep")
    plt.plot(nodes_sweep, s_n(nodes_sweep), 'b--')
    plt.text(28000, 28, '$t_{sweep} = %.6fn^2 + %.6fn + (%.6f)$'%(s[0],s[1],s[2]))
    plt.plot(nodes_brute, times_brute, 'g-', label="Brute")
    plt.plot(nodes_brute, b_n(nodes_brute), 'g--')
    plt.text(17000, 110, '$t_{brute} = %.6fn^2 + %.6fn + (%.6f)$'%(b[0],b[1],b[2]))
    plt.xlabel("Num nodes")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for Average Degree = 16")
    plt.legend(loc=2)
    plt.show()

def plotForVarAvgDeg():
    with open('./report/data/times_cell_deg.pkl', 'r') as f:
        times_cell = pickle.load(f)
    with open('./report/data/avg_degs_cell.pkl', 'r') as f:
        avg_degs_cell = pickle.load(f)

    with open('./report/data/times_sweep_deg.pkl', 'r') as f:
        times_sweep = pickle.load(f)
    with open('./report/data/avg_degs_sweep.pkl', 'r') as f:
        avg_degs_sweep = pickle.load(f)

    c = np.polyfit(avg_degs_cell, times_cell, 1)
    c_n = np.poly1d(c)
    s = np.polyfit(avg_degs_sweep, times_sweep, 1)
    s_n = np.poly1d(s)

    print "cell: t = %.6fn + (%.6f)"%(c[0],c[1])
    print "sweep: t = %.6fn + (%.6f)"%(s[0],s[1])

    plt.plot(avg_degs_cell, times_cell, 'r-', label="Cell")
    plt.plot(avg_degs_cell, c_n(avg_degs_cell), 'r--')
    plt.text(20, 6, '$t_{cell} = %.6fn + (%.6f)$'%(c[0],c[1]))
    plt.plot(avg_degs_sweep, times_sweep, 'b-', label="Sweep")
    plt.plot(avg_degs_sweep, s_n(avg_degs_sweep), 'b--')
    plt.text(12, 56, '$t_{sweep} = %.6fn + (%.6f)$'%(s[0],s[1]))
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
    plt.title("Distribution of Degrees of Nodes, Square")
    plt.show()

def validateIndepSets():
    topology = Square()
    topology.num_nodes = 1000
    topology.avg_deg = 16

    topology.generateNodes()
    topology.findEdges(method="cell")
    topology.colorGraph()

    # check that the color sets are independent
    for i in topology.s_last:
        colors = [topology.node_colors[j] for j in topology.edges[topology.nodes[i]]]
            # colors = [topology.node_colors[topology.s_last.index(j)]]
        if topology.node_colors[i] in colors:
            print "Node shares color with neighbor"

def plotDistributionOfDegreesWhenDel():
    topology = Square()
    topology.num_nodes = 16000
    topology.avg_deg = 32

    topology.generateNodes()
    topology.findEdges(method="cell")
    topology.colorGraph()

    c_orig = Counter([len(v) for v in topology.edges.values()])
    c_del = Counter([topology.deg_when_del[n] for n in topology.nodes])
    labels_orig = [x[0] for x in c_orig.items()]
    labels_del = [x[0] for x in c_del.items()]
    labels = list(set(labels_orig) | set(labels_del))
    orig_deg = [c_orig[l] for l in labels]
    del_deg = [c_del[l] for l in labels]
    indexes = np.arange(len(labels), step=5)
    plt.bar(np.arange(len(labels)), orig_deg, color='b', label="Original Degree")
    plt.bar(np.arange(len(labels)), del_deg, color='r', label="Degree when Deleted in SLVO")
    plt.xticks(indexes + 0.1, labels[::5])
    plt.xlabel("Degree")
    plt.ylabel("Number of Occurances")
    plt.title("Distribution of Degrees of Nodes, Square")
    plt.legend()
    plt.show()

def plotDistributionOfColors():
    topology = Square()
    topology.num_nodes = 16000
    topology.avg_deg = 32

    topology.generateNodes()
    topology.findEdges(method="cell")
    topology.colorGraph()

    c = Counter(topology.node_colors)
    labels = [x[0] for x in c.items()]
    values = [x[1] for x in c.items()]
    indexes = np.arange(len(labels))
    plt.bar(indexes, values)
    plt.xticks(indexes + 0.5, indexes)
    plt.xlabel("Color")
    plt.ylabel("Frequency")
    plt.title("Distribution of Colors, Square")
    plt.show()

def runBenchmarks():
    with open("./report/data/benchmark-data.csv", "w+") as f:
        f.write("Benchmark #,Num. Nodes,Avg. Degree,Topology,Max Degree,Max Deg Deleted,# Color Sets,Max Color Set Size,Terminal Clique Size,Run Time (s)\n")
        n = 0
        tops = {
            'Square': (Square, SQUARE_BENCHMARKS),
            'Disk': (Disk, DISK_BENCHMARKS),
            'Sphere': (Sphere, SPHERE_BENCHMARKS)
        }
        for t in ['Square', 'Disk', 'Sphere']:
            for i in range(len(tops[t][1])):
                n += 1
                topology = tops[t][0]()
                topology.prepBenchmark(i)
                print "nodes: {} | deg: {}".format(topology.num_nodes, topology.avg_deg)
                run_time = time.clock()
                topology.generateNodes()
                topology.findEdges(method="cell")
                topology.colorGraph()
                run_time = time.clock() - run_time
                color_cnt = Counter(topology.node_colors)
                f.write("{},{},{},{},{},".format(n, tops[t][1][i][0], tops[t][1][i][1], t, topology.getMaxDegree()))
                f.write("{},{},{},{},".format(max(topology.deg_when_del.values()), len(set(topology.node_colors)), color_cnt.most_common(1)[0][1], topology.term_clique_size))
                f.write("{0:.3f}\n".format(run_time))
                f.flush()

def main():
    # runForVarNodes()
    # runForVarAvgDeg()
    # plotForVarNodes()
    # plotForVarAvgDeg()
    # plotDistributionOfDegrees()
    # plotDistributionOfDegreesWhenDel()
    # validateIndepSets()
    # plotDistributionOfColors()
    runBenchmarks()

main()
