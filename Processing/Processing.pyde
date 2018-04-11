import random
import time
import math
from collections import Counter
from objects.topology import Square, Disk, Sphere

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 80
AVG_DEG = 10

MAX_NODES_TO_DRAW_EDGES = 8000

RUN_BENCHMARK = False

def setup():
    size(CANVAS_WIDTH, CANVAS_HEIGHT, P3D)
    background(0)

def draw():
    global curr_vis
    
    if curr_vis == 0:
        topology.drawGraph(MAX_NODES_TO_DRAW_EDGES)
    elif curr_vis == 1:
        topology.drawSlvo()
    elif curr_vis == 2:
        topology.drawColoring()

def keyPressed():
    global curr_vis
    global step_size
    
    if key == ' ':
        toggleLooping()
    elif key == 'l':
        incrementVis()
        topology.mightResetCurrNode()
    elif key == 'h':
        decrementVis()
        topology.mightResetCurrNode()
    elif key == 'k':
        topology.incrementCurrNode(step_size)
    elif key == 'j':
        topology.decrementCurrNode(step_size)
    elif key == 'y':
        saveFrame("../report/images/{}-{}.png".format(
                  "slvo" if curr_vis == 1 else "color", topology.curr_node))
    elif key >= '0' and key <= '9':
        step_size = 2**int(key)
        print "New step size:", step_size
    elif key == 'm':
        print "\n---- Help Menu ----"
        print "Use 'hjkl' to move between visualizations"
        print "Press space to pause rotation of the sphere"
        print "Press 'y' to take a screenshot of the current frame"
        print "Entering a number n between 0 and 9 will set the step size to 2^n nodes"

def toggleLooping():
    global is_looping
    if is_looping:
        noLoop()
        is_looping = False
    else:
        loop()
        is_looping = True

def incrementVis():
    global curr_vis
    if curr_vis < 2:
        curr_vis += 1
    background(0)

def decrementVis():
    global curr_vis
    if curr_vis > 0:
        curr_vis -= 1
    background(0)

def main():
    global is_looping
    global curr_vis
    global step_size
    is_looping = True
    curr_vis = 0
    step_size = 1
    
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
    print "Node r: {0:.3f}".format(topology.node_r)
    print "Terminal clique size: {}".format(topology.term_clique_size)
    print "Number of colors: {}".format(len(set(topology.node_colors)))
    print "Max degree when deleted: {}".format(max(topology.deg_when_del.values()))
    color_cnt = Counter(topology.node_colors)
    print "Max color set size: {}  color: {}".format(color_cnt.most_common(1)[0][1],
                                                     color_cnt.most_common(1)[0][0])
    
    run_time = time.clock() - run_time
    print "Run time: {0:.3f} s".format(run_time)

main()
