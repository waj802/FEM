import pandas as pd
import numpy as np
from math import sqrt

def dn_dx_dn_dy_t(dn_dx, dn_dy, size):
    dn_dx_result = []
    dn_dy_result = []
    dn_dx = pd.DataFrame(dn_dx).T
    dn_dy = pd.DataFrame(dn_dy).T

    for x in range(size):
        row_x = pd.DataFrame(dn_dx[x])
        row_x_T = pd.DataFrame(dn_dx[x]).T
        dn_dx_result.append(row_x.dot(row_x_T))

        row_y = pd.DataFrame(dn_dy[x])
        row_y_T = pd.DataFrame(dn_dy[x]).T
        dn_dy_result.append(row_y.dot(row_y_T))

    return dn_dx_result, dn_dy_result


def matrix_dot_det(dn_dx_h, dn_dy_h, det, size):
    for x in range(size):
        array_x = np.array(dn_dx_h[x].values)
        array_x = array_x * float(det[x])
        dn_dx_h[x] = pd.DataFrame(array_x)

        array_y = np.array(dn_dy_h[x].values)
        array_y = array_y * float(det[x])
        dn_dy_h[x] = pd.DataFrame(array_y)

    return dn_dx_h, dn_dy_h


def matrix_dot_k(dn_dx_h, dn_dy_h, k, size):
    def multiply(matrix_x, matrix_y):
        matrix = [[], [], [], []]
        for z in range(4):
            for y in range(4):
                matrix[z].append(round(k * (matrix_x[z][y] + matrix_y[z][y]), 10))
        return matrix

    k_matrix = []
    for i in range(size):
        k_matrix.append(pd.DataFrame(multiply(dn_dx_h[i], dn_dy_h[i])))

    return k_matrix


def final_h_matrix(part_matrix_h, size):
    matrix_h = np.zeros((4, 4))

    if size == 4:
        for i in range(4):
            for j in range(4):
                matrix_h[i][j] = part_matrix_h[0][i][j] + \
                                 part_matrix_h[1][i][j] + \
                                 part_matrix_h[2][i][j] + \
                                 part_matrix_h[3][i][j]
    else:
        for i in range(4):
            for j in range(4):
                matrix_h[i][j] = part_matrix_h[0][i][j] * (5 / 9) * (5 / 9) + \
                                 part_matrix_h[1][i][j] * (8 / 9) * (5 / 9) + \
                                 part_matrix_h[2][i][j] * (5 / 9) * (5 / 9) + \
                                 part_matrix_h[3][i][j] * (5 / 9) * (8 / 9) + \
                                 part_matrix_h[4][i][j] * (8 / 9) * (8 / 9) + \
                                 part_matrix_h[5][i][j] * (5 / 9) * (8 / 9) + \
                                 part_matrix_h[6][i][j] * (5 / 9) * (5 / 9) + \
                                 part_matrix_h[7][i][j] * (8 / 9) * (5 / 9) + \
                                 part_matrix_h[8][i][j] * (5 / 9) * (5 / 9)

    return matrix_h

