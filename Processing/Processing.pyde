import random
import sys
import time
import math
from collections import Counter
from objects.topology import Square, Disk, Sphere

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 20
AVG_DEG = 10

MAX_NODES_TO_DRAW_EDGES = 8000

RUN_BENCHMARK = False

def setup():
    size(CANVAS_WIDTH, CANVAS_HEIGHT, P3D)
    background(0)

def draw():
    global curr_vis
    global draw_domination
    
    if curr_vis == 0:
        topology.drawGraph(MAX_NODES_TO_DRAW_EDGES)
    elif curr_vis == 1:
        topology.drawSlvo()
    elif curr_vis == 2:
        topology.drawColoring()
    elif curr_vis == 3:
        topology.drawPairs(0)
    elif curr_vis == 4:
        topology.drawPairs(1)
    elif curr_vis == 5:
        topology.drawPairs(2)
    elif curr_vis == 6:
        topology.drawPairs(3)
    elif curr_vis == 7:
        topology.drawBackbones(draw_domination)

def keyPressed():
    global curr_vis
    global step_size
    global vis_names
    
    if key == ' ':
        toggleLooping()
    elif key == 'c':
        if curr_vis == 7:
            toggleDrawDomination()
    elif key == 'i':
        topology.switchFgBg()
    elif key == 'l':
        incrementVis()
        topology.mightResetCurrNode()
        print vis_names[curr_vis]
    elif key == 'h':
        decrementVis()
        topology.mightResetCurrNode()
        print vis_names[curr_vis]
    elif key == 'k':
        if curr_vis > 2 and curr_vis < 7:
            topology.incrementCurrPair()
        elif curr_vis == 7:
            topology.incrementCurrBackbone()
        else:
            topology.incrementCurrNode(step_size)
    elif key == 'j':
        if curr_vis > 2 and curr_vis < 7:
            topology.decrementCurrPair()
        elif curr_vis == 7:
            topology.decrementCurrBackbone()
        else:
            topology.decrementCurrNode(step_size)
    elif key == 'y':
        saveFrame("../report/images/{}-#####.png".format(vis_names[curr_vis]))
    elif key >= '0' and key <= '9':
        step_size = 2**int(key)
        print "New step size:", step_size
    elif key == ']':
        step_size = 2*step_size
        print "New step size:", step_size
    elif key == '[':
        step_size = step_size/2
        print "New step size:", step_size
    elif key == 'm':
        print "\n---- Help Menu ----"
        print "Use 'hjkl' to move between visualizations"
        print "Press 'i' to invert the color scheme"
        print "Press 'y' to take a screenshot of the current frame"
        print "Press 'c' to show the coverage of the backbone"
        print "Entering a number n between 0 and 9 will set the step size to 2^n nodes"
        print "Using ']' will double the step size, '[' will half it"
        print "Press space to pause rotation of the sphere"

# def mouseDragged():
#     global topology
#     topology.updateRotation(mouseX, mouseY)

def toggleLooping():
    global is_looping
    if is_looping:
        noLoop()
        is_looping = False
    else:
        loop()
        is_looping = True

def toggleDrawDomination():
    global draw_domination
    if draw_domination:
        draw_domination = False
    else:
        draw_domination = True

def incrementVis():
    global curr_vis
    global topology
    if curr_vis < 7:
        curr_vis += 1
    background(topology.color_bg)

def decrementVis():
    global curr_vis
    global topology
    if curr_vis > 0:
        curr_vis -= 1
    background(topology.color_bg)

def main():
    # sys.setrecursionlimit(32000)
    
    global is_looping
    global draw_domination
    global curr_vis
    global step_size
    global vis_names
    is_looping = True
    draw_domination = False
    curr_vis = 0
    step_size = 1
    vis_names = ["rgg", "slvo", "color", "bipartite", "no-tails", "major-comp", "no-bridge", "backbone"]
    
    global topology
    topology = Square()
    # topology = Disk()
    # topology = Sphere()
    
    topology.num_nodes = NUM_NODES
    topology.avg_deg = AVG_DEG
    topology.canvas_height = CANVAS_HEIGHT
    topology.canvas_width = CANVAS_WIDTH
    
    if RUN_BENCHMARK:
        n_benchmark = 1
        topology.prepBenchmark(n_benchmark)
    
    run_time = time.clock()
    
    topology.generateNodes()
    topology.findEdges(method="cell")
    topology.colorGraph()
    topology.generateBackbones()
    
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
    
    print "\nPress 'm' for the menu"

main()
