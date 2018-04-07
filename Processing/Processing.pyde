import random
import time
import math
from collections import Counter
from objects.topology import Square, Disk, Sphere

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 100
AVG_DEG = 10

MAX_NODES_TO_DRAW_EDGES = 8000

RUN_BENCHMARK = False

def setup():
    size(CANVAS_WIDTH, CANVAS_HEIGHT, P3D)
    background(0)

def draw():
    topology.drawGraph(MAX_NODES_TO_DRAW_EDGES)

def main():
    global topology
    # topology = Square()
    # topology = Disk()
    topology = Sphere()
    
    topology.num_nodes = NUM_NODES
    topology.avg_deg = AVG_DEG
    topology.canvas_height = CANVAS_HEIGHT
    topology.canvas_width = CANVAS_WIDTH
    
    if RUN_BENCHMARK:
        n_benchmark = 0
        topology.prepBenchmark(n_benchmark)
    
    run_time = time.clock()
    
    topology.generateNodes()
    topology.findEdges(method="cell")
    topology.colorGraph()
    
    print "Average degree: {}".format(topology.findAvgDegree())
    print "Min degree: {}".format(topology.getMinDegree())
    print "Max degree: {}".format(topology.getMaxDegree())
    print "Num edges: {}".format(topology.findNumEdges())
    print "Terminal clique size: {}".format(topology.term_clique_size)
    print "Number of colors: {}".format(len(set(topology.node_colors)))
    print "Max degree when deleted: {}".format(max(topology.deg_when_del.values()))
    color_cnt = Counter(topology.node_colors)
    print "Max color set size: {}  color: {}".format(color_cnt.most_common(1)[0][1],
                                                     color_cnt.most_common(1)[0][0])
    
    run_time = time.clock() - run_time
    print "Run time: {0:.3f} s".format(run_time)

main()
