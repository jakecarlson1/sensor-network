import random
import time
from objects.topology import Square, Disk, Sphere

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 1000
AVG_DEG = 16

MAX_NODES_TO_DRAW_EDGES = 8000

def setup():
    size(CANVAS_WIDTH, CANVAS_HEIGHT, P3D)
    background(0)

def draw():
    topology.drawNodes()
    if NUM_NODES < MAX_NODES_TO_DRAW_EDGES:
        topology.drawEdges()
    else:
        topology.drawMinMaxDegNodes()

def main():
    global topology
    # topology = Square()
    # topology = Disk()
    topology = Sphere()
    
    topology.num_nodes = NUM_NODES
    topology.avg_deg = AVG_DEG
    topology.canvas_height = CANVAS_HEIGHT
    topology.canvas_width = CANVAS_WIDTH
    
    run_time = time.clock()
    
    topology.generateNodes()
    topology.findEdges(method="cell")
    
    print "Average degree: {}".format(topology.findAvgDegree())
    print "Min degree: {}".format(topology.getMinDegree())
    print "Max degree: {}".format(topology.getMaxDegree())
    
    run_time = time.clock() - run_time
    print "Run time: {0:.3f} s".format(run_time)

main()