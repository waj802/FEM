import pandas as pd
from math import sqrt
import numpy as np
from element import *
from grid import Grid
from node import *
from matrix_H import *


def dn_dksi_calc(value):
    return np.asarray([-1 * 0.25 * (1 - value),
            1 * 0.25 * (1 - value),
            1 * 0.25 * (1 + value),
            -1 * 0.25 * (1 + value)])

def dn_deta_calc(value):
    return np.asarray([-1 * 0.25 * (1 - value),
            -1 * 0.25 * (1 + value),
            1 * 0.25 * (1 + value),
            1 * 0.25 * (1 - value)])

def shape_function(gauss_value):

    if len(gauss_value) == 2:
        ksi = gauss_value[0]
        eta = gauss_value[1]
        ksi_value = [dn_dksi_calc(-ksi),
                     dn_dksi_calc(-ksi),
                     dn_dksi_calc(ksi),
                     dn_dksi_calc(ksi)]
        eta_value = [dn_deta_calc(-eta),
                     dn_deta_calc(eta),
                     dn_deta_calc(eta),
                     dn_deta_calc(-eta)]
    elif len(gauss_value) == 3:
        ksi = gauss_value[0]
        zero = gauss_value[1]
        eta = gauss_value[2]
        
        ksi_value = [dn_dksi_calc(eta),
                     dn_dksi_calc(eta),
                     dn_dksi_calc(eta),
                     dn_dksi_calc(zero),
                     dn_dksi_calc(zero),
                     dn_dksi_calc(zero),
                     dn_dksi_calc(ksi),
                     dn_dksi_calc(ksi),
                     dn_dksi_calc(ksi)]
        eta_value = [dn_deta_calc(eta),
                     dn_deta_calc(zero),
                     dn_deta_calc(ksi),
                     dn_deta_calc(ksi),
                     dn_deta_calc(zero),
                     dn_deta_calc(eta),
                     dn_deta_calc(eta),
                     dn_deta_calc(zero),
                     dn_deta_calc(ksi)]

    return ksi_value, eta_value


def values_calc(x_y_values, ksi_eta_values):
    return (ksi_eta_values[0] * x_y_values[0]) + (ksi_eta_values[1] * x_y_values[1]) + (ksi_eta_values[2] * x_y_values[2]) + (ksi_eta_values[3] * x_y_values[3])

def jacobian_func(x, y, ksi_values, eta_values):

    if len(ksi_values) == 4:
        jacobian = [[], [], [], []]

        for i in range(4):
            jacobian[i].append(values_calc(x, ksi_values[i]))
            jacobian[i].append(values_calc(y, ksi_values[i]))
            jacobian[i].append(values_calc(x, eta_values[i]))
            jacobian[i].append(values_calc(y, eta_values[i]))
    else:
        jacobian = [[], [], [], [],[], [], [], [], []]

        for i in range(9):
            jacobian[i].append(values_calc(x, ksi_values[i]))
            jacobian[i].append(values_calc(y, ksi_values[i]))
            jacobian[i].append(values_calc(x, eta_values[i]))
            jacobian[i].append(values_calc(y, eta_values[i]))

    jacobian = pd.DataFrame(jacobian)
    jacobian = jacobian.T

    return jacobian

def get_det_of_j(jacobian):
    def calc(jacobian):
        return jacobian[3] * jacobian[0] - jacobian[2] * jacobian[1]

    det = []

    if jacobian.size/4 == 4.0:
        for x in range(4):
            det.append(calc(jacobian[x]))
    else:
        for x in range(9):
            det.append(calc(jacobian[x]))
    det = pd.DataFrame(det)
    det = det.T
    return det

def multiply_matrix_1_det(jacobian, det):
    jacobian = jacobian.T
    for _, column in jacobian.iterrows():
        for x in range(4):
            column[x] = (1 / det[x]) * column[x]

    jacobian = jacobian.T
    return jacobian

def dn_dx_dn_dy(d_n_d_ksi, d_n_d_eta, jacobian):
    for index, col in jacobian.iterrows():
        size = len(col)
    if size == 4:
        dn_dx = [[], [], [], []]
        dn_dy = [[], [], [], []]

        for j in range(4):
            for i in range(4):
                dn_dx[j].append(jacobian[j][0] * d_n_d_ksi[j][i] + jacobian[j][1] * d_n_d_eta[j][i])

            for i in range(4):
                dn_dy[j].append(jacobian[j][2] * d_n_d_ksi[j][i] + jacobian[j][3] * d_n_d_eta[j][i])
    else:
        dn_dx = [[], [], [], [], [], [], [], [], []]
        dn_dy = [[], [], [], [], [], [], [], [], []]

        for j in range(9):
            for i in range(4):
                dn_dx[j].append(jacobian[j][0] * d_n_d_ksi[j][i] + jacobian[j][1] * d_n_d_eta[j][i])

            for i in range(4):
                dn_dy[j].append(jacobian[j][2] * d_n_d_ksi[j][i] + jacobian[j][3] * d_n_d_eta[j][i])


        # for i in range(4):
        #     dn_dx[0].append(jacobian[0][0] * d_n_d_ksi[0][i] + jacobian[0][1] * d_n_d_eta[0][i])
        #     dn_dx[1].append(jacobian[1][0] * d_n_d_ksi[1][i] + jacobian[1][1] * d_n_d_eta[1][i])
        #     dn_dx[2].append(jacobian[2][0] * d_n_d_ksi[2][i] + jacobian[2][1] * d_n_d_eta[2][i])
        #     dn_dx[3].append(jacobian[3][0] * d_n_d_ksi[3][i] + jacobian[3][1] * d_n_d_eta[3][i])
        #     dn_dx[4].append(jacobian[4][0] * d_n_d_ksi[4][i] + jacobian[4][1] * d_n_d_eta[4][i])
        #     dn_dx[5].append(jacobian[5][0] * d_n_d_ksi[5][i] + jacobian[5][1] * d_n_d_eta[5][i])
        #     dn_dx[6].append(jacobian[6][0] * d_n_d_ksi[6][i] + jacobian[6][1] * d_n_d_eta[6][i])
        #     dn_dx[7].append(jacobian[7][0] * d_n_d_ksi[7][i] + jacobian[7][1] * d_n_d_eta[7][i])
        #     dn_dx[8].append(jacobian[8][0] * d_n_d_ksi[8][i] + jacobian[8][1] * d_n_d_eta[8][i])
        #
        # for i in range(4):
        #     dn_dy[0].append(jacobian[0][2] * d_n_d_ksi[0][i] + jacobian[0][3] * d_n_d_eta[0][i])
        #     dn_dy[1].append(jacobian[1][2] * d_n_d_ksi[1][i] + jacobian[1][3] * d_n_d_eta[1][i])
        #     dn_dy[2].append(jacobian[2][2] * d_n_d_ksi[2][i] + jacobian[2][3] * d_n_d_eta[2][i])
        #     dn_dy[3].append(jacobian[3][2] * d_n_d_ksi[3][i] + jacobian[3][3] * d_n_d_eta[3][i])
        #     dn_dy[4].append(jacobian[4][2] * d_n_d_ksi[4][i] + jacobian[4][3] * d_n_d_eta[4][i])
        #     dn_dy[5].append(jacobian[5][2] * d_n_d_ksi[5][i] + jacobian[5][3] * d_n_d_eta[5][i])
        #     dn_dy[6].append(jacobian[6][2] * d_n_d_ksi[6][i] + jacobian[6][3] * d_n_d_eta[6][i])
        #     dn_dy[7].append(jacobian[7][2] * d_n_d_ksi[7][i] + jacobian[7][3] * d_n_d_eta[7][i])
        #     dn_dy[8].append(jacobian[8][2] * d_n_d_ksi[8][i] + jacobian[8][3] * d_n_d_eta[8][i])

        # # ***********************ALT****************************
        # for i in range(4):
        #     # dn_dx[i].append(jacobian[i][0] * d_n_d_ksi[i][j] + jacobian[i][1] * d_n_d_eta[i][j])
        #     dn_dx[0].append(jacobian[0][0] * d_n_d_ksi[0][i] * 5/9 + jacobian[0][1] * d_n_d_eta[0][i] * 5/9)
        #     dn_dx[1].append(jacobian[1][0] * d_n_d_ksi[1][i] * 5/9 + jacobian[1][1] * d_n_d_eta[1][i] * 8/9)
        #     dn_dx[2].append(jacobian[2][0] * d_n_d_ksi[2][i] * 5/9 + jacobian[2][1] * d_n_d_eta[2][i] * 5/9)
        #     dn_dx[3].append(jacobian[3][0] * d_n_d_ksi[3][i] * 8/9 + jacobian[3][1] * d_n_d_eta[3][i] * 5/9)
        #     dn_dx[4].append(jacobian[4][0] * d_n_d_ksi[4][i] * 8/9 + jacobian[4][1] * d_n_d_eta[4][i] * 8/9)
        #     dn_dx[5].append(jacobian[5][0] * d_n_d_ksi[5][i] * 8/9 + jacobian[5][1] * d_n_d_eta[5][i] * 5/9)
        #     dn_dx[6].append(jacobian[6][0] * d_n_d_ksi[6][i] * 5/9 + jacobian[6][1] * d_n_d_eta[6][i] * 5/9)
        #     dn_dx[7].append(jacobian[7][0] * d_n_d_ksi[7][i] * 5/9 + jacobian[7][1] * d_n_d_eta[7][i] * 8/9)
        #     dn_dx[8].append(jacobian[8][0] * d_n_d_ksi[8][i] * 5/9 + jacobian[8][1] * d_n_d_eta[8][i] * 5/9)
        #
        # for i in range(4):
        #     # dn_dy[i].append(jacobian[i][2] * d_n_d_ksi[i][j] + jacobian[i][3] * d_n_d_eta[i][j])
        #     dn_dy[0].append(jacobian[0][2] * d_n_d_ksi[0][i] * 5/9 + jacobian[0][3] * d_n_d_eta[0][i] * 5/9)
        #     dn_dy[1].append(jacobian[1][2] * d_n_d_ksi[1][i] * 5/9 + jacobian[1][3] * d_n_d_eta[1][i] * 8/9)
        #     dn_dy[2].append(jacobian[2][2] * d_n_d_ksi[2][i] * 5/9 + jacobian[2][3] * d_n_d_eta[2][i] * 5/9)
        #     dn_dy[3].append(jacobian[3][2] * d_n_d_ksi[3][i] * 8/9 + jacobian[3][3] * d_n_d_eta[3][i] * 5/9)
        #     dn_dy[4].append(jacobian[4][2] * d_n_d_ksi[4][i] * 8/9 + jacobian[4][3] * d_n_d_eta[4][i] * 8/9)
        #     dn_dy[5].append(jacobian[5][2] * d_n_d_ksi[5][i] * 8/9 + jacobian[5][3] * d_n_d_eta[5][i] * 5/9)
        #     dn_dy[6].append(jacobian[6][2] * d_n_d_ksi[6][i] * 5/9 + jacobian[6][3] * d_n_d_eta[6][i] * 5/9)
        #     dn_dy[7].append(jacobian[7][2] * d_n_d_ksi[7][i] * 5/9 + jacobian[7][3] * d_n_d_eta[7][i] * 8/9)
        #     dn_dy[8].append(jacobian[8][2] * d_n_d_ksi[8][i] * 5/9 + jacobian[8][3] * d_n_d_eta[8][i] * 5/9)

    return dn_dx, dn_dy
