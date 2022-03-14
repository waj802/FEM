import pandas as pd
import numpy as np


def calculate_vector_p(element, alpha, values, temp, size):
    def calc(ksi, eta):
        shape_functions = [0.25 * (1 - ksi) * (1 - eta),
                           0.25 * (1 + ksi) * (1 - eta),
                           0.25 * (1 + ksi) * (1 + eta),
                           0.25 * (1 - ksi) * (1 + eta)]

        return shape_functions

    def calc_value(val):
        array = np.array(val.T)
        array = array * alpha * temp
        return array


    def sum_dot_det(pc_list, val1, val2):
        if size == 4:
            sum = pc_list[0] + pc_list[1]
        else:
            sum = pc_list[0] + pc_list[1] + pc_list[2]

        if abs(val1.x - val2.x) == 0:
            det = abs(val1.y - val2.y)/2
        else:
            det = abs(val1.x - val2.x)/2

        sum = sum * det
        return sum

    def vector_p_3_points(element, values):
        vector_p = np.zeros(4)

        if element.nodes[0].flag == True and element.nodes[1].flag == True:

            first_point = [pd.DataFrame(calc(values[2], - 1)), pd.DataFrame(calc(values[1], - 1)), pd.DataFrame(calc(values[0], -1))]

            first_point[0] = calc_value(first_point[0]) * 5/9
            first_point[1] = calc_value(first_point[1]) * 8/9
            first_point[2] = calc_value(first_point[2]) * 5/9

            vector_p = vector_p + sum_dot_det(np.asarray(first_point), element.nodes[0], element.nodes[1])

        if element.nodes[1].flag == True and element.nodes[2].flag == True:

            second_point = [pd.DataFrame(calc(1, values[2])), pd.DataFrame(calc(1, values[1])), pd.DataFrame(calc(1, values[0]))]

            second_point[0] = calc_value(second_point[0]) * 5/9
            second_point[1] = calc_value(second_point[1]) * 8/9
            second_point[2] = calc_value(second_point[2]) * 5/9

            vector_p = vector_p + sum_dot_det(np.asarray(second_point), element.nodes[1], element.nodes[2])


        if element.nodes[2].flag == True and element.nodes[3].flag == True:

            third_point = [pd.DataFrame(calc(values[0], 1)), pd.DataFrame(calc(values[1], 1)), pd.DataFrame(calc(values[2], 1))]

            third_point[0] = calc_value(third_point[0]) * 5/9
            third_point[1] = calc_value(third_point[1]) * 8/9
            third_point[2] = calc_value(third_point[2]) * 5/9

            vector_p = vector_p + sum_dot_det(np.asarray(third_point), element.nodes[2], element.nodes[3])


        if element.nodes[3].flag == True and element.nodes[0].flag == True:

            fourth_point = [pd.DataFrame(calc(- 1, values[2])), pd.DataFrame(calc(- 1, values[1])), pd.DataFrame(calc(- 1, values[0]))]

            fourth_point[0] = calc_value(fourth_point[0]) * 5/9
            fourth_point[1] = calc_value(fourth_point[1]) * 8/9
            fourth_point[2] = calc_value(fourth_point[2]) * 5/9

            vector_p = vector_p + sum_dot_det(np.asarray(fourth_point), element.nodes[3], element.nodes[0])

        return vector_p
    def vector_p_2_points(element, values):
        vector_p = np.zeros(4)

        if element.nodes[0].flag == True and element.nodes[1].flag == True:

            first_point = [pd.DataFrame(calc(values[1], - 1)), pd.DataFrame(calc(values[0], -1))]

            first_point[0] = calc_value(first_point[0])
            first_point[1] = calc_value(first_point[1])

            vector_p = vector_p + sum_dot_det(first_point, element.nodes[0], element.nodes[1])

        if element.nodes[1].flag == True and element.nodes[2].flag == True:

            second_point = [pd.DataFrame(calc(1, values[1])), pd.DataFrame(calc(1, values[0]))]

            second_point[0] = calc_value(second_point[0])
            second_point[1] = calc_value(second_point[1])

            vector_p = vector_p + sum_dot_det(second_point, element.nodes[1], element.nodes[2])


        if element.nodes[2].flag == True and element.nodes[3].flag == True:

            third_point = [pd.DataFrame(calc(values[0], 1)), pd.DataFrame(calc(values[1], 1))]

            third_point[0] = calc_value(third_point[0])
            third_point[1] = calc_value(third_point[1])

            vector_p = vector_p + sum_dot_det(third_point, element.nodes[2], element.nodes[3])


        if element.nodes[3].flag == True and element.nodes[0].flag == True:

            fourth_point = [pd.DataFrame(calc(- 1, values[0])), pd.DataFrame(calc(- 1, values[1]))]

            fourth_point[0] = calc_value(fourth_point[0])
            fourth_point[1] = calc_value(fourth_point[1])

            vector_p = vector_p + sum_dot_det(fourth_point, element.nodes[3], element.nodes[0])

        return vector_p

    if size == 4:
        vector_p = vector_p_2_points(element, values)
    else:
        vector_p = vector_p_3_points(element, values)

    return vector_p
