import random

CANVAS_HEIGHT = 720
CANVAS_WIDTH = 720

NUM_NODES = 10000

def setup():
    size(CANVAS_WIDTH, CANVAS_HEIGHT, P3D)
    background(0)

def draw():
    for n in range(NUM_NODES):
        ellipse(nodes[n][0]*CANVAS_WIDTH, nodes[n][1]*CANVAS_HEIGHT, 5, 5) 

def generateNodes():
    nodes = []
    for i in range(NUM_NODES):
        nodes.append((random.uniform(0,1), random.uniform(0,1)))

    return nodes

def main():
    global nodes
    nodes = generateNodes()
    
main()