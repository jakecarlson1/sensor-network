import random
import time

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 1000
AVG_DEG = 16

MAX_NODES_TO_DRAW_EDGES = 100

"""
Topology - super class for the shape of the random geometric graph 
"""
class Topology(object):
    
    def __init__(self):
        self.nodes = []
        self.edges = {}
        self.node_r = 0.0
    
    def generateNodes(self):
        print "Method for generating nodes not subclassed"
            
    def findEdges(self, method="brute"):
        self.getRadiusForAverageDegree()
        self.addNodesAsEdgeKeys()
        
        if method == "brute":
            self.bruteForceFindEdges()
    
    def bruteForceFindEdges(self):
        for n in self.nodes:
            for m in self.nodes:
                if n != m and self.distance(n, m) <= self.node_r:
                    self.edges[n].append(m)
    
    def getRadiusForAverageDegree(self):
        print "Method for finding necessary radius for average degree not subclassed"
    
    def addNodesAsEdgeKeys(self):
        self.edges = dict((n, []) for n in self.nodes)
    
    def distance(self, n, m):
        return sqrt((n[0] - m[0])**2+(n[1] - m[1])**2)
    
    def findAvgDegree(self):
        sigma_degree = 0
        for k in self.edges.keys():
            sigma_degree += len(self.edges[k])
        
        return sigma_degree/len(self.edges.keys()) 
    
    def drawNodes(self):
        strokeWeight(2)
        stroke(255)
        fill(255)
        
        for n in range(NUM_NODES):
            ellipse(self.nodes[n][0]*CANVAS_WIDTH, self.nodes[n][1]*CANVAS_HEIGHT, 5, 5)
    
    def drawEdges(self):
        strokeWeight(1)
        stroke(245)
        fill(255)
        
        for n in self.edges.keys():
            for m in self.edges[n]:
                line(n[0]*CANVAS_WIDTH, n[1]*CANVAS_HEIGHT, m[0]*CANVAS_WIDTH, m[1]*CANVAS_HEIGHT)
    
    def drawMinMaxDegNodes(self):
        minDeg = self.edges.keys()[0]
        maxDeg = self.edges.keys()[0]
        
        for k in self.edges.keys():
            if len(self.edges[k]) < len(self.edges[minDeg]):
                minDeg = k
            if len(self.edges[k]) > len(self.edges[maxDeg]):
                maxDeg = k
        
        strokeWeight(1)
        stroke(0,255,0)
        fill(255)
        for n in self.edges[minDeg]:
            line(minDeg[0]*CANVAS_WIDTH, minDeg[1]*CANVAS_HEIGHT, n[0]*CANVAS_WIDTH, n[1]*CANVAS_HEIGHT)
        
        stroke(0,0,255)
        for n in self.edges[maxDeg]:
            line(maxDeg[0]*CANVAS_WIDTH, maxDeg[1]*CANVAS_HEIGHT, n[0]*CANVAS_WIDTH, n[1]*CANVAS_HEIGHT)

"""
Square - inherits from Topology, overloads generateNodes and getRadiusForAverageDegree
for a unit square topology
"""
class Square(Topology):
    
    def __init__(self):
        super(Square, self).__init__()
    
    def generateNodes(self):
        for i in range(NUM_NODES):
            self.nodes.append((random.uniform(0,1), random.uniform(0,1)))
    
    def getRadiusForAverageDegree(self):
        self.node_r = sqrt(AVG_DEG/(NUM_NODES * PI))

"""
Disk - inherits from Topology, overloads generateNodes and getRadiusForAverageDegree
for a unit circle topology
"""
class Disk(Topology):
    
    def __init__(self):
        super(Disk, self).__init__()
    
    def generateNodes(self):
        for i in range(NUM_NODES):
            p = (random.uniform(0,1), random.uniform(0,1))
            while self.distance(p, (0.5,0.5)) > 0.5:
                p = (random.uniform(0,1), random.uniform(0,1))
            self.nodes.append(p)
    
    def getRadiusForAverageDegree(self):
        self.node_r = sqrt((AVG_DEG + 0.0)/NUM_NODES)/2

"""
Sphere - inherits from Topology, overloads generateNodes, getRadiusForAverageDegree,
and distance for a unit sphere topology
"""
class Sphere(Topology):
    
    def __init__(self):
        super(Sphere, self).__init__()
        
    def generateNodes(self):
        for i in range(NUM_NODES):
            p = (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))
            while self.distance(p, (0.5,0.5,0.5)) > 0.5:
                p = (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))
            self.nodes.append(p)

    def getRadiusForAverageDegree(self):
            self.node_r = sqrt((AVG_DEG + 0.0)/NUM_NODES)

    def distance(self, n, m):
        return sqrt((n[0] - m[0])**2+(n[1] - m[1])**2+(n[2] - m[2])**2)

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
    
    run_time = time.clock()
    
    topology.generateNodes()
    topology.findEdges()
    
    print "Average degree: {}".format(topology.findAvgDegree())
    
    run_time = time.clock() - run_time
    print "Run time: {0:.3f} s".format(run_time)
    
main()