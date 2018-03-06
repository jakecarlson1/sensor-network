import random
import math

# benchmarks (num_nodes, avg_deg)
SQUARE_BENCHMARKS = [(1000,32), (8000,64), (16000,32), (64000,64), (64000,128),
                     (128000,64), (128000, 128)]
DISK_BENCHMARKS = [(8000,64), (64000,64), (64000,128)]
SPHERE_BENCHMARKS = [(16000,64), (32000,128), (64000,128)]

"""
Topology - super class for the shape of the random geometric graph
"""
class Topology(object):

    num_nodes = 100
    avg_deg = 0
    canvas_height = 720
    canvas_width = 720

    def __init__(self):
        self.nodes = []
        self.edges = {}
        self.node_r = 0.0
        self.minDeg = ()
        self.maxDeg = ()

    # public funciton for generating nodes of the graph, must be subclassed
    def generateNodes(self):
        print "Method for generating nodes not subclassed"

    # public function for finding edges
    def findEdges(self, method="brute"):
        self._getRadiusForAverageDegree()
        self._addNodesAsEdgeKeys()

        if method == "brute":
            self._bruteForceFindEdges()
        elif method == "sweep":
            self._sweepFindEdges()
        elif method == "cell":
            self._cellFindEdges()
        else:
            print "Find edges method not defined: {}".format(method)

        self._findMinAndMaxDegree()

    # brute force edge detection
    def _bruteForceFindEdges(self):
        for n in self.nodes:
            for m in self.nodes:
                if n != m and self._distance(n, m) <= self.node_r:
                    self.edges[n].append(m)

    # sweep edge detection (2D)
    def _sweepFindEdges(self):
        self.nodes.sort(key=lambda x: x[0])

        for i, n in enumerate(self.nodes):
            search_space = []
            for j in range(1,i+1):
                if abs(n[0] - self.nodes[i-j][0]) <= self.node_r:
                    search_space.append(self.nodes[i-j])
                else:
                    break
            for j in range(1,self.num_nodes-i):
                if abs(n[0] - self.nodes[i+j][0]) <= self.node_r:
                    search_space.append(self.nodes[i+j])
                else:
                    break
            for m in search_space:
                if self._distance(n, m) <= self.node_r:
                    self.edges[n].append(m)

    # cell edge detection (2D)
    def _cellFindEdges(self):
        num_cells = int(1/self.node_r) + 1
        cells = []
        for i in range(num_cells):
            cells.append([[] for j in range(num_cells)])

        for n in self.nodes:
            cells[int(n[0]/self.node_r)][int(n[1]/self.node_r)].append(n)

        for i in range(num_cells):
            for j in range(num_cells):
                for n in cells[i][j]:
                    for c in self._findAdjCells(i, j, num_cells):
                        for m in cells[c[0]][c[1]]:
                            if n != m and self._distance(n, m) <= self.node_r:
                                self.edges[n].append(m)

    # cell edge detection helper function (2D)
    def _findAdjCells(self, i, j, n):
        result = []
        xRange = [(i-1)%n, i, (i+1)%n]
        yRange = [(j-1)%n, j, (j+1)%n]
        for x in xRange:
            for y in yRange:
                result.append((x,y))

        return result

    # function for finding the radius needed for the desired average degree
    # must be subclassed
    def _getRadiusForAverageDegree(self):
        print "Method for finding necessary radius for average degree not subclassed"

    # helper function for findEdges, initializes edges dict
    def _addNodesAsEdgeKeys(self):
        self.edges = dict((n, []) for n in self.nodes)

    # claculates the distance between two nodes (2D)
    def _distance(self, n, m):
        return math.sqrt((n[0] - m[0])**2+(n[1] - m[1])**2)

    # public function for finding the average degree of nodes
    def findAvgDegree(self):
        sigma_degree = 0
        for k in self.edges.keys():
            sigma_degree += len(self.edges[k])

        return sigma_degree/len(self.edges.keys())

    # helper funciton for finding nodes with min and max degree
    def _findMinAndMaxDegree(self):
        self.minDeg = self.edges.keys()[0]
        self.maxDeg = self.edges.keys()[0]

        for k in self.edges.keys():
            if len(self.edges[k]) < len(self.edges[self.minDeg]):
                self.minDeg = k
            if len(self.edges[k]) > len(self.edges[self.maxDeg]):
                self.maxDeg = k

    # public function for getting the minimum degree
    def getMinDegree(self):
        return len(self.edges[self.minDeg])

    # public functino for getting the maximum degree
    def getMaxDegree(self):
        return len(self.edges[self.maxDeg])

    # public function for setting up the benchmark to run, must be subclassed
    def prepBenchmark(self, n):
        print "Method for preparing benchmark not subclassed"

    # public function for drawing the graph
    def drawGraph(self, n_limit):
        self._drawNodes()
        if self.num_nodes <= n_limit:
            self._drawEdges()
        else:
            self._drawMinMaxDegNodes()

    # responsible for drawing the nodes in the canvas
    def _drawNodes(self):
        strokeWeight(2)
        stroke(255)
        fill(255)

        for n in range(self.num_nodes):
            ellipse(self.nodes[n][0]*self.canvas_width, self.nodes[n][1]*self.canvas_height, 5, 5)

    # responsible for drawing the edges in the canavas
    def _drawEdges(self):
        strokeWeight(1)
        stroke(245)
        fill(255)

        for n in self.edges.keys():
            for m in self.edges[n]:
                line(n[0]*self.canvas_width, n[1]*self.canvas_height, m[0]*self.canvas_width, m[1]*self.canvas_height)

    # responsible for drawing the edges of the min and max degree nodes
    def _drawMinMaxDegNodes(self):
        strokeWeight(1)
        stroke(0,255,0)
        fill(255)
        for n in self.edges[self.minDeg]:
            line(self.minDeg[0]*self.canvas_width, self.minDeg[1]*self.canvas_height, n[0]*self.canvas_width, n[1]*self.canvas_height)

        stroke(0,0,255)
        for n in self.edges[self.maxDeg]:
            line(self.maxDeg[0]*self.canvas_width, self.maxDeg[1]*self.canvas_height, n[0]*self.canvas_width, n[1]*self.canvas_height)

    # uses smallest last vertex ordering to color the graph
    def colorGraph(self):
        s_last = self._smallestLastVertexOrdering()

    # constructs a degree structure and determines the smallest last vertex ordering
    def _smallestLastVertexOrdering(self):
        deg_list = [[] for i in range(len(self.edges[self.maxDeg])+1)]

        for n in self.nodes:
            deg_list[len(self.edges[n])].append(n)

        smallest_last_ordering = []

        clique_found = False
        j = len(self.nodes)
        while j > 0:
            # get the current smallest bucket
            curr_bucket = 0
            while curr_bucket < len(deg_list) and len(deg_list[curr_bucket]) == 0:
                curr_bucket += 1

            # if all the remaining nodes are connected we have the terminal clique
            if not clique_found and len(deg_list[curr_bucket]) == j:
                clique_found = True
                smallest_last_ordering.extend(deg_list[curr_bucket])
                self.term_clique_size = curr_bucket
                print "Terminal clique size:", self.term_clique_size

            # get node with smallest degree
            v = deg_list[curr_bucket].pop()
            smallest_last_ordering.append(v)

            # decrement position of nodes that shared an edge with v
            for n in self.edges[v]:
                for i in range(len(deg_list)):
                    to_decrement = True
                    try:
                        deg_list[i].remove(n)
                    except:
                        to_decrement = False

                    if to_decrement:
                        deg_list[i-1].append(n)

            j -= 1

        # reverse list since it was built shortest-first
        return smallest_last_ordering[::-1]

"""
Square - inherits from Topology, overloads generateNodes and _getRadiusForAverageDegree
for a unit square topology
"""
class Square(Topology):

    def __init__(self):
        super(Square, self).__init__()

    # places nodes uniformly in a unit square
    def generateNodes(self):
        for i in range(self.num_nodes):
            self.nodes.append((random.uniform(0,1), random.uniform(0,1)))

    # calculates the radius needed for the requested average degree in a unit square
    def _getRadiusForAverageDegree(self):
        self.node_r = math.sqrt(self.avg_deg/(self.num_nodes * math.pi))

    # gets benchmark setting for square
    def prepBenchmark(self, n):
        self.num_nodes = SQUARE_BENCHMARKS[n][0]
        self.avg_deg = SQUARE_BENCHMARKS[n][1]

"""
Disk - inherits from Topology, overloads generateNodes and _getRadiusForAverageDegree
for a unit circle topology
"""
class Disk(Topology):

    def __init__(self):
        super(Disk, self).__init__()

    # places nodes uniformly in a unit square and regenerates the node if it falls
    # outside of the circle
    def generateNodes(self):
        for i in range(self.num_nodes):
            p = (random.uniform(0,1), random.uniform(0,1))
            while self._distance(p, (0.5,0.5)) > 0.5:
                p = (random.uniform(0,1), random.uniform(0,1))
            self.nodes.append(p)

    # calculates the radius needed for the requested average degree in a unit circle
    def _getRadiusForAverageDegree(self):
        self.node_r = math.sqrt((self.avg_deg + 0.0)/self.num_nodes)/2

    # gets benchmark setting for disk
    def prepBenchmark(self, n):
        self.num_nodes = DISK_BENCHMARKS[n][0]
        self.avg_deg = DISK_BENCHMARKS[n][1]

"""
Sphere - inherits from Topology, overloads generateNodes, _getRadiusForAverageDegree,
and _distance for a unit sphere topology. Also updates the drawGraph function for
a 3D canvas
"""
class Sphere(Topology):

    # adds rotation and node limit variables
    def __init__(self):
        super(Sphere, self).__init__()
        self.rot = (0,math.pi/4,0) # this may move to Topology if rotation is given to the 2D shapes
        # used to control _drawNodes functionality
        self.n_limit = 8000

    # places nodes in a unit cube and projects them onto the surface of the sphere
    def generateNodes(self):
        for i in range(self.num_nodes):
            # equations for uniformly distributing nodes on the surface area of
            # a sphere: http://mathworld.wolfram.com/SpherePointPicking.html
            u = random.uniform(-1,1)
            theta = random.uniform(0, 2*math.pi)
            p = (
                math.sqrt(1 - u**2) * math.cos(theta),
                math.sqrt(1 - u**2) * math.sin(theta),
                u
            )
            self.nodes.append(p)

    # overrides cell for 3D topology, uses 3D mesh of buckets
    def _cellFindEdges(self):
        num_cells = int(1/self.node_r) + 1
        cells = []
        for i in range(num_cells):
            cells.append([[[] for k in range(num_cells)] for j in range(num_cells)])

        for n in self.nodes:
            cells[int(n[0]/self.node_r)][int(n[1]/self.node_r)][int(n[2]/self.node_r)].append(n)

        for i in range(num_cells):
            for j in range(num_cells):
                for k in range(num_cells):
                    for n in cells[i][j][k]:
                        for c in self._findAdjCells(i, j, k, num_cells):
                            for m in cells[c[0]][c[1]][c[2]]:
                                if n != m and self._distance(n, m) <= self.node_r:
                                    self.edges[n].append(m)

    # overrides adjacent cell finding for 3x3 surrounding buckets
    def _findAdjCells(self, i, j, k, n):
        result = []
        xRange = [(i-1)%n, i, (i+1)%n]
        yRange = [(j-1)%n, j, (j+1)%n]
        zRange = [(k-1)%n, k, (k+1)%n]
        for x in xRange:
            for y in yRange:
                for z in zRange:
                    result.append((x,y,z))

        return result

    # calculates the radius needed for the requested average degree in a unit sphere
    def _getRadiusForAverageDegree(self):
            self.node_r = math.sqrt((self.avg_deg + 0.0)/self.num_nodes)*2

    # calculates the distance between two nodes (3D)
    def _distance(self, n, m):
        return math.sqrt((n[0] - m[0])**2+(n[1] - m[1])**2+(n[2] - m[2])**2)

    # gets benchmark setting for sphere
    def prepBenchmark(self, n):
        self.num_nodes = SPHERE_BENCHMARKS[n][0]
        self.avg_deg = SPHERE_BENCHMARKS[n][1]

    # public function for drawing graph, updates node limit if necessary
    def drawGraph(self, n_limit):
        self.n_limit = n_limit
        self._drawNodesAndEdges()

    # responsible for drawing nodes and edges in 3D space
    def _drawNodesAndEdges(self):
        # positions camera
        camera(self.canvas_width/2, self.canvas_height/2, self.canvas_width*-2, 0.5,0.5,0, 0,1,0)

        # updates rotation
        self.rot = (self.rot[0], self.rot[1]-math.pi/100, self.rot[2])

        background(0)
        strokeWeight(2)
        stroke(255)
        fill(255)

        for n in range(self.num_nodes):
            pushMatrix()

            # sets new rotation
            rotateZ(self.rot[2])
            rotateY(-1*self.rot[1])

            # sets drawing origin to current node
            translate((self.nodes[n][0])*self.canvas_width, (self.nodes[n][1])*self.canvas_height, (self.nodes[n][2])*self.canvas_width)

            # places ellipse at origin
            ellipse(0, 0, 10, 10)

            # draw all edges
            if self.num_nodes <= self.n_limit:
                for e in self.edges[self.nodes[n]]:
                    # draws line from origin to neighboring node
                    line(0,0,0, (e[0] - self.nodes[n][0])*self.canvas_width, (e[1] - self.nodes[n][1])*self.canvas_height, (e[2] - self.nodes[n][2])*self.canvas_width)
            # draw edges for min degree node
            elif self.nodes[n] == self.minDeg:
                stroke(0,255,0)
                for e in self.edges[self.nodes[n]]:
                    # draws line from origin to neighboring node
                    line(0,0,0, (e[0] - self.nodes[n][0])*self.canvas_width, (e[1] - self.nodes[n][1])*self.canvas_height, (e[2] - self.nodes[n][2])*self.canvas_width)
                stroke(255)
            # draw edges for max degree node
            elif self.nodes[n] == self.maxDeg:
                stroke(0,0,255)
                for e in self.edges[self.nodes[n]]:
                    # draws line from origin to neighboring node
                    line(0,0,0, (e[0] - self.nodes[n][0])*self.canvas_width, (e[1] - self.nodes[n][1])*self.canvas_height, (e[2] - self.nodes[n][2])*self.canvas_width)
                stroke(255)

            popMatrix()
