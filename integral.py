from math import sqrt

class point_structure:
    wp = []
    nodes = []
    def add(self, wpp, ndd):
        self.wp += wpp
        self.nodes += ndd

def funx(x):
    return 5*x*x+3*x+6

def funxy(x, y):
    return 5*x*x*y*y+3*x*y+6

def Gaussx(wp, nodes):
    total = 0
    for p in range(len(wp)):
        x = nodes[p]
        total += wp[p] * funx(x)
    return total

def Gaussxy(wp, nodes):
    total = 0
    for p1 in range(len(wp)):
        for p2 in range(len(wp)):
            x = nodes[p1]
            y = nodes[p2]
            total += wp[p1] * wp[p2] * funxy(x, y)
    return total



class Element4_2D(object):

    def dn_dksi(self, eta):
            return -0.25*(1-eta)

    def dn_deta(self, ksi):
            return -0.25*(1-ksi)

    def calc(self, nodes):
            matrix = []
            for p1 in range(len(nodes)):
                for p2 in range(len(nodes)):
                    matrix[p1][p2] = dn_dksi(nodes[p2])
            return matrix

if __name__ == "__main__":
    p = point_structure()
    p.add([1,1], [-1/sqrt(3), 1/sqrt(3)])
    # p.add([5/9,8/9,5/9], [-sqrt(3/5), 0, sqrt(3/5)])

    # print(Gaussx(p.wp, p.nodes))
    # print(Gaussxy(p.wp, p.nodes))

    e = Element4_2D()
    print(e.calc([-1/sqrt(3), 1/sqrt(3)]))