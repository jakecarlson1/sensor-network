import random
import math
import time
from collections import deque

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
        self.slvo = []
        self.deg_when_del = {}
        self.node_colors = []
        self.pairs = []
        self.no_tails = []
        self.major_comps = []
        self.clean_pairs = []
        self.backbones = []
        self.curr_node = 0
        self.curr_pair = 0
        self.curr_backbone = 0

        self.rot = (0,0,0)
        self.color_bg = 0
        self.color_fg = 255

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
        for i, n in enumerate(self.nodes):
            for j, m in enumerate(self.nodes):
                if i != j and self._distance(n, m) <= self.node_r:
                    self.edges[n].append(j)

    # sweep edge detection
    def _sweepFindEdges(self):
        self.nodes.sort(key=lambda x: x[0])

        for i, n in enumerate(self.nodes):
            search_space = []
            for j in range(1,self.num_nodes-i):
                if abs(n[0] - self.nodes[i+j][0]) <= self.node_r:
                    search_space.append(i+j)
                else:
                    break
            for j in search_space:
                if self._distance(n, self.nodes[j]) <= self.node_r:
                    self.edges[n].append(j)
                    self.edges[self.nodes[j]].append(i)

    # cell edge detection
    def _cellFindEdges(self):
        num_cells = int(1/self.node_r) + 1
        cells = []
        for i in range(num_cells):
            cells.append([[] for j in range(num_cells)])

        for i, n in enumerate(self.nodes):
            cells[int(n[0]/self.node_r)][int(n[1]/self.node_r)].append(i)

        for i in range(num_cells):
            for j in range(num_cells):
                for n_i in cells[i][j]:
                    for c in self._findAdjCells(i, j, num_cells):
                        for m_i in cells[c[0]][c[1]]:
                            if self._distance(self.nodes[n_i], self.nodes[m_i]) <= self.node_r:
                                self.edges[self.nodes[n_i]].append(m_i)
                                self.edges[self.nodes[m_i]].append(n_i)
                    for m_i in cells[i][j]:
                        if self._distance(self.nodes[n_i], self.nodes[m_i]) <= self.node_r and n_i != m_i:
                            self.edges[self.nodes[n_i]].append(m_i)

    # cell edge detection helper function
    def _findAdjCells(self, i, j, n):
        adj_cells = [(1,-1), (0,1), (1,1), (1,0)]
        return (((i+x[0])%n,(j+x[1])%n) for x in adj_cells)

    # function for finding the radius needed for the desired average degree
    # must be subclassed
    def _getRadiusForAverageDegree(self):
        print "Method for finding necessary radius for average degree not subclassed"

    # helper function for findEdges, initializes edges dict
    def _addNodesAsEdgeKeys(self):
        self.edges = {n:[] for n in self.nodes}

    # claculates the distance between two nodes (2D)
    def _distance(self, n, m):
        return math.sqrt((n[0] - m[0])**2+(n[1] - m[1])**2)

    # public function for finding the number of edges
    def findNumEdges(self):
        sigma_edges = 0
        for k in self.edges.keys():
            sigma_edges += len(self.edges[k])

        return sigma_edges/2

    # public function for finding the average degree of nodes
    def findAvgDegree(self):
        return 2*self.findNumEdges()/self.num_nodes

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
        self._drawNodes(self.nodes)
        if self.num_nodes <= n_limit:
            self._drawEdges(self.nodes)
        else:
            self._drawMinMaxDegNodes()

    # responsible for drawing the nodes in the canvas
    def _drawNodes(self, node_list):
        strokeWeight(2)
        stroke(self.color_fg)
        fill(self.color_fg)

        for n in node_list:
            ellipse(n[0]*self.canvas_width, n[1]*self.canvas_height, 5, 5)

    # responsible for drawing the edges in the canavas
    def _drawEdges(self, node_list):
        strokeWeight(1)
        stroke(245)
        fill(self.color_fg)

        s = set(node_list)

        for n in node_list:
            for m_i in self.edges[n]:
                if self.nodes[m_i] in s:
                    line(n[0]*self.canvas_width, n[1]*self.canvas_height, self.nodes[m_i][0]*self.canvas_width, self.nodes[m_i][1]*self.canvas_height)

    # responsible for drawing the edges of the min and max degree nodes
    def _drawMinMaxDegNodes(self):
        strokeWeight(1)
        stroke(0,self.color_fg,0)
        fill(self.color_fg)
        for n_i in self.edges[self.minDeg]:
            line(self.minDeg[0]*self.canvas_width, self.minDeg[1]*self.canvas_height, self.nodes[n_i][0]*self.canvas_width, self.nodes[n_i][1]*self.canvas_height)

        stroke(0,0,self.color_fg)
        for n_i in self.edges[self.maxDeg]:
            line(self.maxDeg[0]*self.canvas_width, self.maxDeg[1]*self.canvas_height, self.nodes[n_i][0]*self.canvas_width, self.nodes[n_i][1]*self.canvas_height)

    # uses smallest last vertex ordering to color the graph
    def colorGraph(self):
        self.slvo, self.deg_when_del = self._smallestLastVertexOrdering()
        self.node_colors = self._assignNodeColors(self.slvo)
        self.color_map = self._mapColorsToRGB(self.node_colors)

    # constructs a degree structure and determines the smallest last vertex ordering
    def _smallestLastVertexOrdering(self):
        deg_sets = {l:set() for l in range(len(self.edges[self.maxDeg])+1)}
        deg_when_del = {n:len(self.edges[n]) for n in self.nodes}

        for i, n in enumerate(self.nodes):
            deg_sets[deg_when_del[n]].add(i)

        smallest_last_ordering = []

        clique_found = False
        j = len(self.nodes)
        while j > 0:
            # get the current smallest bucket
            curr_bucket = 0
            while len(deg_sets[curr_bucket]) == 0:
                curr_bucket += 1

            # if all the remaining nodes are connected we have the terminal clique
            if not clique_found and len(deg_sets[curr_bucket]) == j:
                clique_found = True
                self.term_clique_size = curr_bucket

            # get node with smallest degree
            v_i = deg_sets[curr_bucket].pop()
            smallest_last_ordering.append(v_i)

            # decrement position of nodes that shared an edge with v
            for n_i in (n_i for n_i in self.edges[self.nodes[v_i]] if n_i in deg_sets[deg_when_del[self.nodes[n_i]]]):
                deg_sets[deg_when_del[self.nodes[n_i]]].remove(n_i)
                deg_when_del[self.nodes[n_i]] -= 1
                deg_sets[deg_when_del[self.nodes[n_i]]].add(n_i)

            j -= 1

        # reverse list since it was built shortest-first
        return smallest_last_ordering[::-1], deg_when_del

    # assigns the colors to nodes given in a smallest-last vertex ordering as a parallel array
    def _assignNodeColors(self, slvo):
        colors = [-1 for _ in range(len(slvo))]
        for i in slvo:
            adj_colors = set([colors[j] for j in self.edges[self.nodes[i]]])
            color = 0
            while color in adj_colors:
                color += 1
            colors[i] = color

        return colors

    # generates random color codes for each color set and returns them in a dictionary
    def _mapColorsToRGB(self, color_list):
        s = set(color_list)
        color_map = {}
        while len(s) > 0:
            c = s.pop()
            color_map[c] = (random.randint(0,255), random.randint(0,255), random.randint(0,255))

        return color_map

    # draw nodes as they are removed in smallest-last vertex ordering
    def drawSlvo(self):
        l = [self.nodes[i] for i in self.slvo[0:self.num_nodes - self.curr_node]]
        self._drawNodes(l)
        self._drawEdges(l)

    # increments curr_node, used to limit the number of nodes drawn
    def incrementCurrNode(self, s):
        if self.curr_node + s <= self.num_nodes:
            self.curr_node += s
            background(self.color_bg)
        elif self.curr_node != self.num_nodes:
            self.curr_node = self.num_nodes
            background(self.color_bg)

    # decrements curr_node, used to limit the number of nodes drawn
    def decrementCurrNode(self, s):
        if self.curr_node - s >= 0:
            self.curr_node -= s
            background(self.color_bg)
        elif self.curr_node != 0:
            self.curr_node = 0
            background(self.color_bg)

    # used to reset curr node if all nodes have been drawn and the method changes
    def mightResetCurrNode(self):
        if self.curr_node == self.num_nodes:
            curr_node = 0
            background(self.color_bg)

    # increments curr_backbone, used to draw different backbones
    def incrementCurrPair(self):
        if self.curr_pair < len(self.pairs) - 1:
            self.curr_pair += 1
            background(self.color_bg)

    # decrements curr_backbone, used to draw different backbones
    def decrementCurrPair(self):
        if self.curr_pair > 0:
            self.curr_pair -= 1
            background(self.color_bg)

    # increments curr_backbone, used to draw different backbones
    def incrementCurrBackbone(self):
        if self.curr_backbone < len(self.backbones) - 1:
            self.curr_backbone += 1
            background(self.color_bg)

    # decrements curr_backbone, used to draw different backbones
    def decrementCurrBackbone(self):
        if self.curr_backbone > 0:
            self.curr_backbone -= 1
            background(self.color_bg)

    # switch foreground and background colors
    def switchFgBg(self):
        self.color_fg, self.color_bg = self.color_bg, self.color_fg

    # # update the rotation of the drawing
    # def updateRotation(self, x, y):
    #     # self.rot = (self.rot[0], self.rot[1]-math.pi/100, self.rot[2])
    #     # self.rot = (x*math.cos(self.rot[0])*math.pi/500, self.rot[1], self.rot[2])
    #     self.rot = (self.rot[0], x*math.cos(self.rot[1])*math.pi/1000, self.rot[2])
    #     # rotateX(self.rot[0])
    #     # rotateZ(self.rot[2])
    #     # rotateY(-1*self.rot[1])

    # used to draw the graph with the nodes colored
    def drawColoring(self):
        l = [self.nodes[i] for i in self.slvo[0:self.curr_node]]
        self._drawNodes(l)
        self._applyColors(self.slvo[0:self.curr_node])
        self._drawEdges(l)

    # places colors on the nodes
    def _applyColors(self, node_i_list):
        strokeWeight(5)

        num_colors = max(self.node_colors)

        for n_i in node_i_list:
            c = self.color_map[self.node_colors[n_i]]
            stroke(c[0], c[1], c[2])
            fill(c[0], c[1], c[2])
            ellipse(self.nodes[n_i][0]*self.canvas_width, self.nodes[n_i][1]*self.canvas_height, 5, 5)

    # public function for pairing the independent sets and picking the largest backbones
    def generateBackbones(self):
        # pair four largest independent sets
        self.pairs = self._pairIndependentSets(self.node_colors)

        # delete minor components and tails
        self.no_tails, self.major_comps, self.clean_pairs = self._cleanPairs(self.pairs)

        # pick two backbones of largest size

        # calculate domination

    # pairs the four largest independent color sets
    def _pairIndependentSets(self, color_list):
        # the first four color sets should be the largest (slvo)
        indep_sets = [set() for _ in range(4)]

        for i, n in enumerate(self.nodes):
            if self.node_colors[i] < 4:
                indep_sets[self.node_colors[i]].add(i)

        # return combinations of sets (union)
        return [s1 | s2 for i, s1 in enumerate(indep_sets) for s2 in indep_sets[i+1:]]

    # removes the minor components and tails from the bipartite subgraphs
    def _cleanPairs(self, bipartites):
        no_tails = []
        major_comps = []
        results = []
        for b in bipartites:
            # remove the tails and save the graph for visualization
            b = self._removeTails(b)
            no_tails.append(b)

            # use BFS to get the major component
            major_comp = self._bfs(b)
            major_comps.append(major_comp)

            # use DFS to remove bridges
            backbone = self._removeBridges(major_comp)
            results.append(backbone)


        return no_tails, major_comps, results

    # remove tails from bipartite, very similar to smallest-last vertex ordering
    def _removeTails(self, bipartite):
        bipartite = bipartite.copy()
        # build graph representation
        points = list(bipartite)
        deg_sets = {l:set() for l in range(len(self.edges[self.maxDeg])+1)}
        deg_map = {n_i:len([e_i for e_i in self.edges[self.nodes[n_i]] if e_i in bipartite]) for n_i in points}

        for i, n in enumerate(self.nodes):
            if i in bipartite:
                deg_sets[deg_map[i]].add(i)

        # remove nodes with zero or one edge until there are no tails
        while len(deg_sets[0]) > 0 or len(deg_sets[1]) > 0:
            to_remove = deg_sets[0] | deg_sets[1]
            deg_sets[0] = set()
            deg_sets[1] = set()

            for n_i in list(to_remove):
                for e_i in [e_i for e_i in self.edges[self.nodes[n_i]] if e_i in bipartite]:
                    if e_i in deg_sets[deg_map[e_i]]:
                        deg_sets[deg_map[e_i]].remove(e_i)
                        deg_map[e_i] -= 1
                        deg_sets[deg_map[e_i]].add(e_i)

                bipartite.remove(n_i)

        return bipartite

    # use BFS to find the major component
    def _bfs(self, bipartite, rm_edges=None):
        points = list(bipartite)
        # used to index into the points array
        index_to_local = {n_i:i for i, n_i in enumerate(points)}
        # used to index into the nodes array
        index_to_global = {i:n_i for i, n_i in enumerate(points)}
        visited = [0 for _ in points]
        visits = []
        components = []

        while 0 in visited:
            visit = 1

            queue = deque()
            root = visited.index(0)
            queue.append(root)
            visited[root] = 1
            # builds a set for the points in each component
            components.append(set([index_to_global[root]]))

            while len(queue) > 0:
                curr = queue.pop()

                for e in [index_to_local[e] for e in self.edges[self.nodes[points[curr]]] if e in bipartite]:
                    if rm_edges != None and (e in rm_edges and curr in rm_edges):
                        continue
                    if visited[e] == 0:
                        visit += 1
                        queue.append(e)
                        components[-1].add(index_to_global[e])
                        visited[e] = 1

            visits.append(visit)

        return components[visits.index(max(visits))]

    # removes all bridges and minor blocks from major component
    # algorithm: https://e-maxx-eng.appspot.com/graph/bridge-searching.html
    def _removeBridges(self, major_comp):
        points = list(major_comp)
        # used to index into the points array
        index_to_local = {n_i:i for i, n_i in enumerate(points)}
        # used to index into the nodes array
        index_to_global = {i:n_i for i, n_i in enumerate(points)}
        visited = [0 for _ in points]
        bridge_nodes = set()
        tin = [-1 for _ in points]
        fup = [-1 for _ in points]
        visit = 0

        for i, p in enumerate(points):
            if visited[i] == 0:
                self._dfs(major_comp, points, i, p, index_to_local, visited, bridge_nodes, tin, fup, visit)

        return self._bfs(major_comp, bridge_nodes)

    # use DFS to find bridges
    def _dfs(self, comp, points, i, p, index_to_local, visited, bridge_nodes, tin, fup, visit, to=-1):
        visited[i] = 1
        tin[i] = visit
        fup[i] = visit
        visit += 1
        for e in [index_to_local[e] for e in self.edges[self.nodes[p]] if e in comp]:
            if e == to:
                continue
            if visited[e] == 1:
                fup[i] = min(fup[i], tin[e])
            else:
                self._dfs(comp, points, e, points[e], index_to_local, visited, bridge_nodes, tin, fup, visit, to=i)
                fup[i] = min(fup[i], fup[e])
                if fup[e] > tin[i]:
                    if i not in bridge_nodes:
                        bridge_nodes.add(i)
                    if e not in bridge_nodes:
                        bridge_nodes.add(e)

    # public function for drawing the color set pairs
    def drawPairs(self, mode=0):
        l_i = []
        if mode == 0:
            l_i = list(self.pairs[self.curr_pair])
        elif mode == 1:
            l_i = list(self.no_tails[self.curr_pair])
        elif mode == 2:
            l_i = list(self.major_comp[self.curr_pair])
        elif mode == 3:
            l_i = list(self.clean_pairs[self.curr_pair])
        l_n = [self.nodes[i] for i in l_i]
        self._drawNodes(l_n)
        self._applyColors(l_i)
        self._drawEdges(l_n)

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
        self._drawNodesAndEdges(self.nodes)

    # responsible for drawing nodes and edges in 3D space
    def _drawNodesAndEdges(self, node_list):
        # positions camera
        camera(self.canvas_width/2, self.canvas_height/2, self.canvas_width*-2, 0.5,0.5,0, 0,1,0)

        # updates rotation
        self.rot = (self.rot[0], self.rot[1]-math.pi/100, self.rot[2])

        background(self.color_bg)
        strokeWeight(2)
        stroke(self.color_fg)
        fill(self.color_fg)

        s = set(node_list)

        for n in node_list:
            pushMatrix()

            # sets new rotation
            rotateZ(self.rot[2])
            rotateY(-1*self.rot[1])

            # sets drawing origin to current node
            translate(n[0]*self.canvas_width, n[1]*self.canvas_height, n[2]*self.canvas_width)

            # places ellipse at origin
            ellipse(0, 0, 10, 10)

            # draw all edges
            if self.num_nodes <= self.n_limit:
                for e_i in self.edges[n]:
                    if self.nodes[e_i] in s:
                        e = self.nodes[e_i]
                        # draws line from origin to neighboring node
                        line(0,0,0, (e[0] - n[0])*self.canvas_width, (e[1] - n[1])*self.canvas_height, (e[2] - n[2])*self.canvas_width)
            # draw edges for min degree node
            elif n == self.minDeg:
                stroke(0,self.color_fg,0)
                for e_i in self.edges[n]:
                    e = self.nodes[e_i]
                    # draws line from origin to neighboring node
                    line(0,0,0, (e[0] - n[0])*self.canvas_width, (e[1] - n[1])*self.canvas_height, (e[2] - n[2])*self.canvas_width)
                stroke(self.color_fg)
            # draw edges for max degree node
            elif n == self.maxDeg:
                stroke(0,0,self.color_fg)
                for e_i in self.edges[n]:
                    e = self.nodes[e_i]
                    # draws line from origin to neighboring node
                    line(0,0,0, (e[0] - n[0])*self.canvas_width, (e[1] - n[1])*self.canvas_height, (e[2] - n[2])*self.canvas_width)
                stroke(self.color_fg)

            popMatrix()

    # draw nodes as they are removed in smallest-last vertex ordering
    def drawSlvo(self):
        l = [self.nodes[i] for i in self.slvo[0:self.num_nodes - self.curr_node]]
        self._drawNodesAndEdges(l)

    # used to draw the graph with the nodes colored
    def drawColoring(self):
        l = [self.nodes[i] for i in self.slvo[0:self.curr_node]]
        self._drawNodesAndEdges(l)
        self._applyColors(self.slvo[0:self.curr_node])

    # places colors on the nodes
    def _applyColors(self, node_i_list):
        strokeWeight(2)

        num_colors = max(self.node_colors)

        for n_i in node_i_list:
            c = self.color_map[self.node_colors[n_i]]
            stroke(c[0], c[1], c[2])
            fill(c[0], c[1], c[2])

            pushMatrix()

            # sets new rotation
            rotateZ(self.rot[2])
            rotateY(-1*self.rot[1])

            # sets drawing origin to current node
            translate(self.nodes[n_i][0]*self.canvas_width, self.nodes[n_i][1]*self.canvas_height, self.nodes[n_i][2]*self.canvas_width)

            # places ellipse at origin
            ellipse(0, 0, 10, 10)

            popMatrix()

    # public function for drawing the color set pairs
    def drawPairs(self, mode=0):
        l_i = []
        if mode == 0:
            l_i = list(self.pairs[self.curr_pair])
        elif mode == 1:
            l_i = list(self.no_tails[self.curr_pair])
        elif mode == 2:
            l_i = list(self.major_comps[self.curr_pair])
        elif mode == 3:
            l_i = list(self.clean_pairs[self.curr_pair])
        l_n = [self.nodes[i] for i in l_i]
        self._drawNodesAndEdges(l_n)
        self._applyColors(l_i)
