import random
import time
from objects.topology import Square, Disk, Sphere

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 1000
AVG_DEG = 16

MAX_NODES_TO_DRAW_EDGES = 16000

# """
# Topology - super class for the shape of the random geometric graph 
# """
# class Topology(object):
    
#     def __init__(self):
#         self.nodes = []
#         self.edges = {}
#         self.node_r = 0.0
    
#     def generateNodes(self):
#         print "Method for generating nodes not subclassed"
            
#     def findEdges(self, method="brute"):
#         self.getRadiusForAverageDegree()
#         self.addNodesAsEdgeKeys()
        
#         if method == "brute":
#             self.bruteForceFindEdges()
#         elif method == "sweep":
#             self.sweepFindEdges()
#         elif method == "cell":
#             self.cellFindEdges()
    
#     def bruteForceFindEdges(self):
#         for n in self.nodes:
#             for m in self.nodes:
#                 if n != m and self.distance(n, m) <= self.node_r:
#                     self.edges[n].append(m)
    
#     def sweepFindEdges(self):
#         self.nodes.sort(key=lambda x: x[0])
        
#         for i, n in enumerate(self.nodes):
#             search_space = []
#             for j in range(1,i+1):
#                 if abs(n[0] - self.nodes[i-j][0]) <= self.node_r:
#                     search_space.append(self.nodes[i-j])
#                 else:
#                     break
#             for j in range(1,NUM_NODES-i):
#                 if abs(n[0] - self.nodes[i+j][0]) <= self.node_r:
#                     search_space.append(self.nodes[i+j])
#                 else:
#                     break
#             for m in search_space:
#                 if self.distance(n, m) <= self.node_r:
#                     self.edges[n].append(m)
    
#     def cellFindEdges(self):
#         num_cells = int(1/self.node_r) + 1
#         cells = []
#         for i in range(num_cells):
#             cells.append([[] for j in range(num_cells)])
        
#         for n in self.nodes:
#             cells[int(n[0]/self.node_r)][int(n[1]/self.node_r)].append(n)
        
#         for i in range(num_cells):
#             for j in range(num_cells):
#                 for n in cells[i][j]:
#                     for c in self.findPosAdjCells(i, j, num_cells):
#                         for m in cells[c[0]][c[1]]:
#                             if n != m and self.distance(n, m) <= self.node_r:
#                                 self.edges[n].append(m)

#     def findPosAdjCells(self, i, j, n):
#         result = []
#         xRange = [(i-1)%n, i, (i+1)%n]
#         yRange = [(j-1)%n, j, (j+1)%n]
#         for x in xRange:
#             for y in yRange:
#                 result.append((x,y))
        
#         return result
    
#     def getRadiusForAverageDegree(self):
#         print "Method for finding necessary radius for average degree not subclassed"
    
#     def addNodesAsEdgeKeys(self):
#         self.edges = dict((n, []) for n in self.nodes)
    
#     def distance(self, n, m):
#         return sqrt((n[0] - m[0])**2+(n[1] - m[1])**2)
    
#     def findAvgDegree(self):
#         sigma_degree = 0
#         for k in self.edges.keys():
#             sigma_degree += len(self.edges[k])
        
#         return sigma_degree/len(self.edges.keys()) 
    
#     def drawNodes(self):
#         strokeWeight(2)
#         stroke(255)
#         fill(255)
        
#         for n in range(NUM_NODES):
#             ellipse(self.nodes[n][0]*CANVAS_WIDTH, self.nodes[n][1]*CANVAS_HEIGHT, 5, 5)
    
#     def drawEdges(self):
#         strokeWeight(1)
#         stroke(245)
#         fill(255)
        
#         for n in self.edges.keys():
#             for m in self.edges[n]:
#                 line(n[0]*CANVAS_WIDTH, n[1]*CANVAS_HEIGHT, m[0]*CANVAS_WIDTH, m[1]*CANVAS_HEIGHT)
    
#     def drawMinMaxDegNodes(self):
#         minDeg = self.edges.keys()[0]
#         maxDeg = self.edges.keys()[0]
        
#         for k in self.edges.keys():
#             if len(self.edges[k]) < len(self.edges[minDeg]):
#                 minDeg = k
#             if len(self.edges[k]) > len(self.edges[maxDeg]):
#                 maxDeg = k
        
#         strokeWeight(1)
#         stroke(0,255,0)
#         fill(255)
#         for n in self.edges[minDeg]:
#             line(minDeg[0]*CANVAS_WIDTH, minDeg[1]*CANVAS_HEIGHT, n[0]*CANVAS_WIDTH, n[1]*CANVAS_HEIGHT)
        
#         stroke(0,0,255)
#         for n in self.edges[maxDeg]:
#             line(maxDeg[0]*CANVAS_WIDTH, maxDeg[1]*CANVAS_HEIGHT, n[0]*CANVAS_WIDTH, n[1]*CANVAS_HEIGHT)

# """
# Square - inherits from Topology, overloads generateNodes and getRadiusForAverageDegree
# for a unit square topology
# """
# class Square(Topology):
    
#     def __init__(self):
#         super(Square, self).__init__()
    
#     def generateNodes(self):
#         for i in range(NUM_NODES):
#             self.nodes.append((random.uniform(0,1), random.uniform(0,1)))
    
#     def getRadiusForAverageDegree(self):
#         self.node_r = sqrt(AVG_DEG/(NUM_NODES * PI))

# """
# Disk - inherits from Topology, overloads generateNodes and getRadiusForAverageDegree
# for a unit circle topology
# """
# class Disk(Topology):
    
#     def __init__(self):
#         super(Disk, self).__init__()
    
#     def generateNodes(self):
#         for i in range(NUM_NODES):
#             p = (random.uniform(0,1), random.uniform(0,1))
#             while self.distance(p, (0.5,0.5)) > 0.5:
#                 p = (random.uniform(0,1), random.uniform(0,1))
#             self.nodes.append(p)
    
#     def getRadiusForAverageDegree(self):
#         self.node_r = sqrt((AVG_DEG + 0.0)/NUM_NODES)/2

# """
# Sphere - inherits from Topology, overloads generateNodes, getRadiusForAverageDegree,
# and distance for a unit sphere topology
# """
# class Sphere(Topology):
    
#     def __init__(self):
#         super(Sphere, self).__init__()

#     # TODO: fix this as it fills the volume of the sphere
#     def generateNodes(self):
#         for i in range(NUM_NODES):
#             p = (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))
#             while self.distance(p, (0.5,0.5,0.5)) > 0.5:
#                 p = (random.uniform(0,1), random.uniform(0,1), random.uniform(0,1))
#             self.nodes.append(p)

#     def getRadiusForAverageDegree(self):
#             self.node_r = sqrt((AVG_DEG + 0.0)/NUM_NODES)

#     def distance(self, n, m):
#         return sqrt((n[0] - m[0])**2+(n[1] - m[1])**2+(n[2] - m[2])**2)
    
#     # TODO: overload drawNodes and drawEdges for a 3D object, add rotation

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
    topology = Square()
    # topology = Disk()
    # topology = Sphere()
    
    topology.num_nodes = NUM_NODES
    topology.avg_deg = AVG_DEG
    topology.canvas_height = CANVAS_HEIGHT
    topology.canvas_width = CANVAS_WIDTH
    
    run_time = time.clock()
    
    topology.generateNodes()
    topology.findEdges(method="cell")
    
    print "Average degree: {}".format(topology.findAvgDegree())
    
    run_time = time.clock() - run_time
    print "Run time: {0:.3f} s".format(run_time)

main()