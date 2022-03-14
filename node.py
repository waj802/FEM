from dataclasses import dataclass
import math

@dataclass
class Node(object):

    def __init__(self, x, y, node_id):
        super(Node, self).__init__()
        self.x = x
        self.y = y
        self.id = node_id
        self.flag = False

    def check_node(self, element):
        self.flag = self.set_flag(element.H, element.B)

    def __repr__(self):
        return f'({self.x} : {self.y})'

    def set_flag(self, h, b):

        if math.isclose(self.x, 0, abs_tol=0.001) or math.isclose(self.x, b, abs_tol=0.001):
            return True
        elif math.isclose(self.y, 0, abs_tol=0.001) or math.isclose(self.y, h, abs_tol=0.001):
            return True
        return False