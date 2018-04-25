import sys
import time
import matplotlib.pyplot as plt
import numpy as np
import pickle
from Processing.objects.topology import *
from collections import Counter

def varryNumNodes(top=Square):
    print "-------- {} --------".format(top)
    times = []
    nodes = range(4000,68000,4000)
    for n in nodes:
        topology = top()
        topology.num_nodes = n
        topology.avg_deg = 16

        time_0 = time.clock()

        topology.generateNodes()
        topology.findEdges(method="cell")
        topology.colorGraph()
        topology.generateBackbones()

        time_0 = time.clock() - time_0
        times.append(time_0)

        print "{}: {} s".format(n, time_0)

        if time_0 > 90:
            break

    return times, nodes[:len(times)]

def varryAvgDeg(top=Square):
    print "-------- {} --------".format(top)
    times = []
    avg_degs = range(2,34,2)
    for d in avg_degs:
        topology = top()
        topology.num_nodes = 32000
        topology.avg_deg = d

        time_0 = time.clock()

        topology.generateNodes()
        topology.findEdges(method="cell")
        topology.colorGraph()
        topology.generateBackbones()

        time_0 = time.clock() - time_0
        times.append(time_0)

        print "{}: {} s".format(d, time_0)

        if time_0 > 90:
            break

    return times, avg_degs[:len(times)]

def runForVarNodes():
    times_square, nodes_square = varryNumNodes(top=Square)
    with open('./report/data/times_square.pkl', 'w') as f:
        pickle.dump(times_square, f)
    with open('./report/data/nodes_square.pkl', 'w') as f:
        pickle.dump(nodes_square, f)

    times_disk, nodes_disk = varryNumNodes(top=Disk)
    with open('./report/data/times_disk.pkl', 'w') as f:
        pickle.dump(times_disk, f)
    with open('./report/data/nodes_disk.pkl', 'w') as f:
        pickle.dump(nodes_disk, f)

    times_sphere, nodes_sphere = varryNumNodes(top=Sphere)
    with open('./report/data/times_sphere.pkl', 'w') as f:
        pickle.dump(times_sphere, f)
    with open('./report/data/nodes_sphere.pkl', 'w') as f:
        pickle.dump(nodes_sphere, f)

def runForVarAvgDeg():
    times_square, avg_degs_square = varryAvgDeg(top=Square)
    with open('./report/data/square_times_deg.pkl', 'w') as f:
        pickle.dump(times_square, f)
    with open('./report/data/square_avg_degs.pkl', 'w') as f:
        pickle.dump(avg_degs_square, f)

    times_disk, avg_degs_disk = varryAvgDeg(top=Disk)
    with open('./report/data/disk_times_deg.pkl', 'w') as f:
        pickle.dump(times_disk, f)
    with open('./report/data/disk_avg_degs.pkl', 'w') as f:
        pickle.dump(avg_degs_disk, f)

    times_sphere, avg_degs_sphere = varryAvgDeg(top=Sphere)
    with open('./report/data/sphere_times_deg.pkl', 'w') as f:
        pickle.dump(times_sphere, f)
    with open('./report/data/sphere_avg_degs.pkl', 'w') as f:
        pickle.dump(avg_degs_sphere, f)

def plotForVarNodes():
    with open('./report/data/times_square.pkl', 'r') as f:
        times_square = pickle.load(f)
    with open('./report/data/nodes_square.pkl', 'r') as f:
        nodes_square = pickle.load(f)

    with open('./report/data/times_disk.pkl', 'r') as f:
        times_disk = pickle.load(f)
    with open('./report/data/nodes_disk.pkl', 'r') as f:
        nodes_disk = pickle.load(f)

    with open('./report/data/times_sphere.pkl', 'r') as f:
        times_sphere = pickle.load(f)
    with open('./report/data/nodes_sphere.pkl', 'r') as f:
        nodes_sphere = pickle.load(f)

    sq = np.polyfit(nodes_square, times_square, 1)
    sq_n = np.poly1d(sq)
    d = np.polyfit(nodes_disk, times_disk, 1)
    d_n = np.poly1d(d)
    sp = np.polyfit(nodes_sphere, times_sphere, 1)
    sp_n = np.poly1d(sp)

    print "square: t = %.6fn + (%.6f)"%(sq[0],sq[1])
    print "disk: t = %.6fn + (%.6f)"%(d[0],d[1])
    print "sphere: t = %.6fn + (%.6f)"%(sp[0],sp[1])

    plt.plot(nodes_square, times_square, 'b-', label="Square")
    plt.plot(nodes_square, sq_n(nodes_square), 'b--')
    plt.text(40000, 5, '$t_{square} = %.6fn + (%.6f)$'%(sq[0],sq[1]))
    plt.plot(nodes_disk, times_disk, 'g-', label="Disk")
    plt.plot(nodes_disk, d_n(nodes_disk), 'g--')
    plt.text(40000, 1.5, '$t_{disk} = %.6fn + (%.6f)$'%(d[0],d[1]))
    plt.plot(nodes_sphere, times_sphere, 'r-', label="Sphere")
    plt.plot(nodes_sphere, sp_n(nodes_sphere), 'r--')
    plt.text(14000, 16, '$t_{sphere} = %.6fn + (%.6f)$'%(sp[0],sp[1]))
    plt.xlabel("Num nodes")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for Average Degree = 16")
    plt.legend(loc=2)
    plt.show()

def plotForVarAvgDeg():
    with open('./report/data/square_times_deg.pkl', 'r') as f:
        times_square = pickle.load(f)
    with open('./report/data/square_avg_degs.pkl', 'r') as f:
        avg_degs_square = pickle.load(f)

    with open('./report/data/disk_times_deg.pkl', 'r') as f:
        times_disk = pickle.load(f)
    with open('./report/data/disk_avg_degs.pkl', 'r') as f:
        avg_degs_disk = pickle.load(f)

    with open('./report/data/sphere_times_deg.pkl', 'r') as f:
        times_sphere = pickle.load(f)
    with open('./report/data/sphere_avg_degs.pkl', 'r') as f:
        avg_degs_sphere = pickle.load(f)

    sq = np.polyfit(avg_degs_square, times_square, 1)
    sq_n = np.poly1d(sq)
    d = np.polyfit(avg_degs_disk, times_disk, 1)
    d_n = np.poly1d(d)
    sp = np.polyfit(avg_degs_sphere, times_sphere, 1)
    sp_n = np.poly1d(sp)

    print "square: t = %.6fn + (%.6f)"%(sq[0],sq[1])
    print "disk: t = %.6fn + (%.6f)"%(d[0],d[1])
    print "sphere: t = %.6fn + (%.6f)"%(sp[0],sp[1])

    plt.plot(avg_degs_square, times_square, 'b-', label="Square")
    plt.plot(avg_degs_square, sq_n(avg_degs_square), 'b--')
    plt.text(18, 1, '$t_{square} = %.6fn + (%.6f)$'%(sq[0],sq[1]))
    plt.plot(avg_degs_disk, times_disk, 'g-', label="Disk")
    plt.plot(avg_degs_disk, d_n(avg_degs_disk), 'g--')
    plt.text(18, 4, '$t_{disk} = %.6fn + (%.6f)$'%(d[0],d[1]))
    plt.plot(avg_degs_sphere, times_sphere, 'r-', label="Sphere")
    plt.plot(avg_degs_sphere, sp_n(avg_degs_sphere), 'r--')
    plt.text(8, 15, '$t_{sphere} = %.6fn + (%.6f)$'%(sp[0],sp[1]))
    plt.xlabel("Average Degree")
    plt.ylabel("Run time (s)")
    plt.title("Run Time for |V| = 32,000")
    plt.legend(loc=2)
    plt.show()

def plotDistributionOfDegrees(topology, top=None, sep=5):
    c = Counter([len(v) for v in topology.edges.values()])
    labels = [x[0] for x in c.items()]
    values = [x[1] for x in c.items()]
    indexes = np.arange(len(labels))
    plt.bar(indexes, values)
    plt.xticks(indexes[::sep] + 0.5, labels[::sep])
    plt.xlabel("Degree")
    plt.ylabel("Number of Occurances")
    plt.title("Distribution of Degrees of Nodes{}".format(", {}, |V| = {}, A = {}".format(top, topology.num_nodes, topology.avg_deg) if top != None else ""))
    plt.show()

def validateIndepSets(topology):
    # check that the color sets are independent
    valid = True
    for i in topology.slvo:
        colors = [topology.node_colors[j] for j in topology.edges[topology.nodes[i]]]
        if topology.node_colors[i] in colors:
            valid = False

    print "Valid coloring" if valid else "Node shares color with neighbor"

def plotDistributionOfDegreesWhenDel(topology, top=None, sep=5):
    indexes = np.arange(len(topology.nodes), step=sep)
    # verticies, degrees
    plt.plot(np.arange(len(topology.nodes)), [len(topology.edges[topology.nodes[n_i]]) for n_i in topology.slvo], 'g', label="Original")
    plt.plot(np.arange(len(topology.nodes)), [topology.deg_when_del[topology.nodes[n_i]] for n_i in topology.slvo], 'b', label="Deleted")
    plt.xlabel("Vertex")
    plt.ylabel("Degree")
    plt.title("Vertex Degree Plots in SLVO Ordering{}".format(", {}, |V| = {}, A = {}".format(top, topology.num_nodes, topology.avg_deg) if top != None else ""))
    plt.legend()
    plt.show()

def plotDistributionOfColors(topology, top=None, sep=5):
    c = Counter(topology.node_colors)
    labels = [x[0] for x in c.items()]
    values = [x[1] for x in c.items()]
    indexes = np.arange(len(labels))
    plt.bar(indexes, values)
    plt.xticks(indexes[::sep] + 0.5, indexes[::sep])
    plt.xlabel("Color")
    plt.ylabel("Frequency")
    plt.title("Distribution of Colors{}".format(", {}, |V| = {}, A = {}".format(top, topology.num_nodes, topology.avg_deg) if top != None else ""))
    plt.show()

def runBenchmarks(graphs=False):
    with open("./report/data/benchmark-data-1.csv", "w+") as f1:
        with open("./report/data/benchmark-data-2.csv", "w+") as f2:
            with open("./report/data/benchmark-data-3.csv", "w+") as f3:
                f1.write("Benchmark,Order,A,Topology,r,Size,Realized A,Max Deg,Min Deg,Run Time (s)\n")
                f2.write("Benchmark,Max Deg Deleted,Color Sets,Largest Color Set,Terminal Clique Size\n")
                f3.write("Benchmark,B1 Colors,B1 Order,B1 Size,B1 Domination,B1 Faces,B2 Colors,B2 Order,B2 Size,B2 Domination,B2 Faces\n")
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
                        topology.generateBackbones()

                        run_time = time.clock() - run_time

                        color_cnt = Counter(topology.node_colors)
                        b1_colors = list(set([str(topology.node_colors[j]) for j in list(topology.backbones[0])]))
                        b2_colors = list(set([str(topology.node_colors[j]) for j in list(topology.backbones[1])]))

                        f1.write("{},{},{},{},".format(n, tops[t][1][i][0], tops[t][1][i][1], t))
                        f1.write("{0:.3f},".format(topology.node_r))
                        f1.write("{},{},{},{},".format(topology.findNumEdges(), topology.findAvgDegree(), topology.getMaxDegree(), topology.getMinDegree()))
                        f1.write("{0:.3f}\n".format(run_time))
                        f2.write("{},{},{},{},{}\n".format(n, max(topology.deg_when_del.values()), len(set(topology.node_colors)), color_cnt.most_common(1)[0][1], topology.term_clique_size))
                        f3.write("{},{},{},{},".format(n, " \& ".join(b1_colors), topology.backbones_meta[0][0], topology.backbones_meta[0][1]))
                        f3.write("{0:.3f},".format(topology.backbones_meta[0][2]))
                        f3.write("{},{},{},{},".format(topology.num_faces[0] if t == "Sphere" else "", " \& ".join(b2_colors), topology.backbones_meta[1][0], topology.backbones_meta[1][1]))
                        f3.write("{0:.3f},".format(topology.backbones_meta[1][2]))
                        f3.write("{}\n".format(topology.num_faces[1] if t == "Sphere" else ""))
                        f1.flush()
                        f2.flush()
                        f3.flush()

                        validateIndepSets(topology)

                        if graphs:
                            sep = 5
                            if topology.avg_deg >= 128:
                                sep = 25
                            plotDistributionOfDegrees(topology, t, sep)
                            plotDistributionOfDegreesWhenDel(topology, t, sep)
                            plotDistributionOfColors(topology, t)

def main():
    # topology = Square()
    # topology.num_nodes = 32000
    # topology.avg_deg = 16
    #
    # topology.generateNodes()
    # topology.findEdges(method="cell")
    # topology.colorGraph()
    # topology.generateBackbones()

    # runForVarNodes()
    # runForVarAvgDeg()
    # plotForVarNodes()
    # plotForVarAvgDeg()
    # plotDistributionOfDegrees(topology)
    # plotDistributionOfDegreesWhenDel(topology, sep=1)
    # validateIndepSets(topology)
    # plotDistributionOfColors(topology, sep=1)
    # sys.setrecursionlimit(8000)
    # runBenchmarks()

    tops = {
        'Square': (Square, SQUARE_BENCHMARKS),
        'Disk': (Disk, DISK_BENCHMARKS),
        'Sphere': (Sphere, SPHERE_BENCHMARKS)
    }
    for t in ['Square', 'Disk', 'Sphere']:
        for i in range(len(tops[t][1])):
            with open('./report/data/{}_{}.pkl'.format(t, i), 'r') as f:
                topology = pickle.load(f)

                sep = 5
                if topology.avg_deg >= 128:
                    sep = 25
                # plotDistributionOfDegrees(topology, t, sep)
                # plotDistributionOfDegreesWhenDel(topology, t, sep)
                # plotDistributionOfColors(topology, t)

main()
