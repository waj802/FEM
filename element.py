from dataclasses import dataclass
from matrix_H import *

@dataclass
class Element(object):
    def __init__(self, ID, H, B, nH, nB, grid):
        super(Element, self).__init__()
        self.ID = ID
        self.H = H
        self.B = B
        self.nH = nH
        self.nB = nB
        self.grid = grid
        self.setPrecision = 10
        self.x = -1
        self.y = -1
        self.deltaY = -1
        self.deltaX = -1
        self.nodes = []
        self.det = 0
        self.getCoordinates()
        self.generateNodes()
        self.check_status_of_flag()

    def __str__(self):
        return f'Element ID: {self.ID} Nodes: {self.nodes}'


    def getCoordinates(self):  #ustalanie koordynatow na podstawie zmiany x i y oraz id elementu
        self.deltaX = float(str(self.B / (self.nB - 1))[:self.setPrecision])
        self.deltaY = float(str(self.H / (self.nH - 1))[:self.setPrecision])

        loop = self.ID
        deltaX = self.deltaX
        deltaY = self.deltaY
        self.x = deltaX / 2
        self.y = deltaY / 2
        while True:
            if deltaX * loop > self.B:
                loop -= self.nB - 1
                self.y += deltaY
            else:
                self.x = deltaX * loop - deltaX / 2
                break
        self.x = round(self.x, 10)
        self.y = round(self.y, 10)




    def generateNodes(self):   #ustalanie wezlow dla elementu oraz walidacja w getNode
        self.nodes = [self.grid.getNode(self.x - (self.deltaX / 2),
                                        self.y - (self.deltaY / 2)),
                      self.grid.getNode(self.x + (self.deltaX / 2),
                                        self.y - (self.deltaY / 2)),
                      self.grid.getNode(self.x + (self.deltaX / 2),
                                        self.y + (self.deltaY / 2)),
                      self.grid.getNode(self.x - (self.deltaX / 2),
                                        self.y + (self.deltaY / 2))]

    def check_status_of_flag(self):
        for n in self.nodes:
            n.check_node(self)

