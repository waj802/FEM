from dataclasses import dataclass
from element import Element
from node import Node
import numpy as np
import math

@dataclass
class Grid(object):
    def __init__(self, H, B, nH, nB):
        super(Grid, self).__init__()
        self.H = H
        self.B = B
        self.nH = nH
        self.nB = nB
        self.elements = []
        self.nodes: Node = []

        self.global_h_matrix = np.zeros((nH * nB, nB * nB))
        self.global_c_matrix = np.zeros((nH * nB, nB * nB))
        self.global_p_vector = np.zeros((nH * nB))

        self.generateNodes(nH, nB)
        self.generateElements(nH, nB)


    def generateNodes(self, nH, nB):   #ustalanie koordynatow na podstawie zmiany x i y oraz id wezla
        deltaX = self.B / (self.nB - 1)
        deltaY = self.H / (self.nH - 1)

        for i in range(nH * nB):
            x: float = 0
            y: float = 0
            loop = i
            while True:
                if deltaX * loop <= self.B:
                    x += deltaX * loop
                    break
                else:
                    loop -= self.nB
                    y += deltaY
            self.nodes.append(Node(x, y, i))


    def generateElements(self, nH, nB):
        nE = (nH - 1) * (nB - 1)    #max liczba elementow
        elementsArray = []
        for e in range(nE):
            elementsArray.append(Element(e + 1, self.H, self.B, nH, nB, self))
        self.elements = elementsArray

    def getNode(self, x0, y0):  #sprawdzanie koordynatow x i y oraz niwelowanie bledu przyblizenia
        for node in self.nodes:
            if math.isclose(node.x, x0, abs_tol=0.001) and math.isclose(node.y, y0, abs_tol=0.001):
            # if abs(node.x - x0) <= 0.0001 and abs(node.y - y0) <= 0.0001:
                return node

    def elementsPrint(self):
        for e in self.elements:
            print(e)