import numpy as np
import pandas as pd

from grid import Grid
from matrix_H import *
from jacobian import *
from main import *
from matrix_hbc import *
from vector_P import *
from matrix_C import *


def run_matrix_h_grid(element, values, k):
    x = [element.nodes[0].x, element.nodes[1].x, element.nodes[2].x, element.nodes[3].x]
    y = [element.nodes[0].y, element.nodes[1].y, element.nodes[2].y, element.nodes[3].y]


    dn_dksi, dn_deta = shape_function(values)

    jacobian = jacobian_func(x, y, dn_dksi, dn_deta)

    element.det = get_det_of_j(jacobian)

    jacobian = multiply_matrix_1_det(jacobian, element.det)

    dn_dx, dn_dy = dn_dx_dn_dy(dn_dksi, dn_deta, jacobian)

    for _, col in jacobian.iterrows():
        size = len(col)

    dn_dx_h, dn_dy_h = dn_dx_dn_dy_t(dn_dx, dn_dy, size)

    dn_dx_h, dn_dy_h = matrix_dot_det(dn_dx_h, dn_dy_h, element.det, size)

    h_matrix = final_h_matrix(matrix_dot_k(dn_dx_h, dn_dy_h, k, size), size)

    return h_matrix, size

def calc_matrix_h_global(element, matrix_h_lok, grid):
    nodes = []
    for n in element.nodes:
        nodes.append(n.id)

    for i in range(4):
        for j in range(4):
            grid.global_h_matrix[nodes[i], nodes[j]] += matrix_h_lok[i,j]


def calc_matrix_c_global(element, matrix_c_lok, grid):
    nodes = []
    for n in element.nodes:
        nodes.append(n.id)

    for i in range(4):
        for j in range(4):
            grid.global_c_matrix[nodes[i], nodes[j]] += matrix_c_lok[i,j]


def calc_vector_p_global(element, vector_p_lok, grid):
    nodes = []
    for n in element.nodes:
        nodes.append(n.id)

    try:
        for val, idx in zip(vector_p_lok[0], range(grid.nB*grid.nH)):
            grid.global_p_vector[nodes[idx]] += val
    except:
        for val, idx in zip(vector_p_lok, range(grid.nB*grid.nH)):
            grid.global_p_vector[nodes[idx]] += val



def run(H, B, nH, nB, t0_temp, time_step, time_final, k, alfa, tot, density, spec_heat, *args):

    #schemat 2 punktowy
    ksi = 1/sqrt(3)
    eta = -1/sqrt(3)
    values = [ksi, eta]


    #schemat 3 punktowy
    # ksi = sqrt(3/5)
    # zero = 0
    # eta = -sqrt(3/5)
    # values = [ksi,zero,eta]

    g = Grid(H,B,nH,nB)

    t0 = np.zeros(g.nB * g.nH)
    for i in range(g.nB * g.nH):
        t0[i] = t0_temp

    for i in range(int(time_final/time_step)):
        for element in g.elements:
            matrix_h, size = run_matrix_h_grid(element, values, k)
            matrix_hbc = calculate_matrix_hbc(element, alfa, values, size)
            matrix_h_lok = matrix_hbc+matrix_h
            calc_matrix_h_global(element, matrix_h_lok, g)

            vectorp = calculate_vector_p(element, alfa, values, tot, size)
            calc_vector_p_global(element, vectorp, g)

            matrix_c_lok = np.array(calculate_matrix_c(element.det, values, spec_heat, density, size))
            calc_matrix_c_global(element, matrix_c_lok, g)
        matrix_c_divide = g.global_c_matrix/time_step
        h_result = g.global_h_matrix+matrix_c_divide

        print("Macierz H: ")
        print(pd.DataFrame(g.global_h_matrix))
        print("Macierz C: ")
        print(pd.DataFrame(g.global_c_matrix))
        print("Wektor P: ")
        print(pd.DataFrame(g.global_p_vector))
        print("Wektor final: ")
        print(pd.DataFrame(h_result))
        # exit()

        p_result = g.global_p_vector + (matrix_c_divide.dot(t0))
        print(pd.DataFrame(p_result))
        exit()

        t1 = np.linalg.solve(h_result, p_result)

        print(f'Step time: {i*time_step+time_step} Temp min: {round(np.min(t1),3)}  Temp max: {round(np.max(t1),3)}')

        t0 = t1
        


if __name__ == "__main__":

    H = 0.1  #wysokosc
    B = 0.1  #szerokosc
    nH = 4   #liczba wezlow po H
    nB = 4   #liczba wezlow po B
    # nH = 31   #liczba wezlow po H
    # nB = 31   #liczba wezlow po B

    # load("test.txt")

    run(H=H, B=B, nH=nH, nB=nB, t0_temp=100, time_step=50, time_final=500, k=25, alfa=300, tot=1200, density=7800, spec_heat=700)
    # run(H=H, B=B, nH=nH, nB=nB, t0_temp=100, time_step=1, time_final=100, k=25, alfa=300, tot=1200, density=7800, spec_heat=700)

    # run(H=H, B=B, nH=nH, nB=nB, t0_temp=100, time_step=50, time_final=500)
    # run(H=H, B=B, nH=nH, nB=nB, t0_temp=100, time_step=1, time_final=100)