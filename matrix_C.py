import pandas as pd
import numpy as np


def calculate_matrix_c(det, values, c, p, size):

    def calc_N(values):
        if size == 4:
            ksi = values[0]
            eta = values[1]
        else:
            ksi = values[0]
            zero = values[1]
            eta = values[2]
        shape_function_values = []

        def calc_ksi_eta(ksi, eta):
            shape_values_local = []
            shape_values_local.append(0.25 * (1 - ksi) * (1 - eta))
            shape_values_local.append(0.25 * (1 + ksi) * (1 - eta))
            shape_values_local.append(0.25 * (1 + ksi) * (1 + eta))
            shape_values_local.append(0.25 * (1 - ksi) * (1 + eta))

            return np.asarray(shape_values_local)

        if size == 4:
            shape_function_values.append(calc_ksi_eta(eta, eta))
            shape_function_values.append(calc_ksi_eta(ksi, eta))
            shape_function_values.append(calc_ksi_eta(ksi, ksi))
            shape_function_values.append(calc_ksi_eta(eta, ksi))
        else:
            shape_function_values.append(calc_ksi_eta(eta, eta))
            shape_function_values.append(calc_ksi_eta(zero, eta))
            shape_function_values.append(calc_ksi_eta(ksi, eta))
            shape_function_values.append(calc_ksi_eta(eta, zero))
            shape_function_values.append(calc_ksi_eta(zero, zero))
            shape_function_values.append(calc_ksi_eta(ksi, zero))
            shape_function_values.append(calc_ksi_eta(eta, ksi))
            shape_function_values.append(calc_ksi_eta(zero, ksi))
            shape_function_values.append(calc_ksi_eta(ksi, ksi))        

        return pd.DataFrame(shape_function_values)

    shape_function = calc_N(values)
    shape_function = shape_function.T

    matrix_sum = []

    for i in range(size):
        x = pd.DataFrame(shape_function[i])
        x_trans = pd.DataFrame(shape_function[i]).T

        result = x.dot(x_trans)
        result = np.array(result)
        result = p * c * result * float(det[i])
        matrix_sum.append(result)

    matrix_c = np.zeros((4, 4))

    if size == 4:
        for i in range(4):
            for j in range(4):
                matrix_c[i][j] = matrix_sum[0][i][j] + \
                                 matrix_sum[1][i][j] + \
                                 matrix_sum[2][i][j] + \
                                 matrix_sum[3][i][j]
    else:
        for i in range(4):
            for j in range(4):
                matrix_c[i][j] = matrix_sum[0][i][j] * (5 / 9) * (5 / 9) + \
                                 matrix_sum[1][i][j] * (8 / 9) * (5 / 9) + \
                                 matrix_sum[2][i][j] * (5 / 9) * (5 / 9) + \
                                 matrix_sum[3][i][j] * (5 / 9) * (8 / 9) + \
                                 matrix_sum[4][i][j] * (8 / 9) * (8 / 9) + \
                                 matrix_sum[5][i][j] * (5 / 9) * (8 / 9) + \
                                 matrix_sum[6][i][j] * (5 / 9) * (5 / 9) + \
                                 matrix_sum[7][i][j] * (8 / 9) * (5 / 9) + \
                                 matrix_sum[8][i][j] * (5 / 9) * (5 / 9)

    return pd.DataFrame(matrix_c)