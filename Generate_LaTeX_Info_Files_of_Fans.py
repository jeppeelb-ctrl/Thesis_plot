import re
import numpy as np
import sys
import pyperclip
from io import StringIO
import math
import string
# This block is to generate the LaTeX code for the data of a fan
# It generates: TikZ code for the fan
#               Displays the generators and their representation
#               Computes and displays the anticanonical divisor
#               Rescales the anticanoncial divisor
#               Checks and displays if the anticanonical divisor is Kähler
#               Computes and displays the intersection matrix of the fan

# Redirect stdout to capture print statements
output_buffer = StringIO()
sys.stdout = output_buffer


def define_cones(list_of_rays):
    '''
    This function creates a list of the cones of the fan

    list_of_rays: is a list of the rays of the fan

    Returns: a list of the cones of the fan
    '''

    cones_of_fan = []

    for i in range(len(list_of_rays)):
        if i == len(list_of_rays) - 1:
            cones_of_fan.append(Cone([list_of_rays[i], list_of_rays[0]]))
        else:
            cones_of_fan.append(Cone([list_of_rays[i], list_of_rays[i + 1]]))

    return cones_of_fan


def create_fan(rays_of_fan, cones_of_fan):
    '''
    This function creates the fan given the rays and cones of the fan

    rays_of_fan: is a list of the rays of the fan of the form fx. (0,1)

    cones_of_fan: is a list of the cones of the fan

    Returns: the fan with rays and cones corresponding to rays_of_fan and cones_of_fan respectively
    '''

    fan = Fan(cones=cones_of_fan, rays=rays_of_fan)

    return fan


def print_everything(input_strings):

    def sort_from_x(rays):
        # Split the string into a list of vectors
        vectors = rays.split(';')

        # Convert each vector to a tuple of integers
        points = []
        for vec in vectors:
            x, y = map(int, vec.split(','))
            points.append((x, y))

        # Sort the points anticlockwise starting from (1, 0)
        start_point = (1, 0)
        sorted_points = sorted(points, key=lambda point: math.atan2(point[1], point[0]))

        # Rotate the list so that (1, 0) is first
        start_index = sorted_points.index(start_point)
        sorted_points = sorted_points[start_index:] + sorted_points[:start_index]


        # Convert the sorted points back to the string format
        fan_str_rays = []
        for ray in sorted_points:
            ray = list(ray)
            fan_str_rays.append(ray)
        return fan_str_rays


    def get_dr(point):
        dr = sqrt(point[0]**2 + point[1]**2)
        return dr


    def get_dx(point):
        dx = point[0]
        return dx


    def get_dy(point):
        dy = point[1]
        return dy


    def get_x(point):
        x1 = (abs(get_dx(point)) * sqrt(get_dr(point)**2)) / (get_dr(point) ** 2)
        return x1


    def get_y(point):
        y1 = (abs(get_dy(point)) * sqrt(get_dr(point)**2)) / (get_dr(point) ** 2)
        return y1

    def get_coordinate(point):
        if point[0] < 0:
            x = -round(get_x(point), 4)
        elif point[0] >= 0:
            x = round(get_x(point), 4)
        if point[1] < 0:
            y = -round(get_y(point), 4)
        elif point[1] >= 0:
            y = round(get_y(point), 4)
        return [x, y]


    def get_angle(point1, point2):
        dot_product = sum(i*j for i, j in zip(point1, point2))
        norm_point1 = sqrt(sum(i**2 for i in point1))
        norm_point2 = sqrt(sum(i**2 for i in point2))
        cos_theta = dot_product / (norm_point1 * norm_point2)
        angle_rad = math.acos(cos_theta)
        angle_deg = math.degrees(angle_rad)
        return angle_deg


    def get_cummulative_angles(rays):
        cumm_angles = []
        for i in range(len(rays)):
            if i == 0:
                ang = get_angle(rays[i], rays[i + 1])
                cumm_angles.append([0, ang])
            elif i == len(rays) - 1:
                prev_ang = cumm_angles[i - 1][1]
                new_ang = 360
                cumm_angles.append([prev_ang, new_ang])
            else:
                ang = get_angle(rays[i], rays[i + 1])
                prev_ang = cumm_angles[i - 1][1]
                new_ang = prev_ang + ang
                cumm_angles.append([prev_ang, new_ang])
        return cumm_angles


    def get_cones(rays):
        cones = []
        for i in range(len(rays)):
            if i == len(rays) - 1:
                cones.append([rays[i], rays[0]])
            else:
                cones.append([rays[i], rays[i + 1]])
        return cones


    def get_mid_line(cone):
        x = (cone[0][0] + cone[1][0]) / 2
        y = (cone[0][1] + cone[1][1]) / 2
        return [x,y]


    def get_angle_rot(vector):
        x, y = vector
        return math.atan2(y, x)

    def convert_ray_list(rays):
        start_vector = [0, 1]

        sorted_vectors = sorted(rays, key=lambda vector: -get_angle_rot(vector))

        sorted_vectors = sorted_vectors[sorted_vectors.index(start_vector):] + sorted_vectors[:sorted_vectors.index(start_vector)]

        return sorted_vectors


    def extract_coefficients(s):
        # Split the string into terms
        terms = re.findall(r'([+-]?\s*\d*[*]?z\d+|[+-]?\s*z\d+)', s)
        result = {}
        for term in terms:
            # Remove any whitespace
            term = term.replace(' ', '')
            # Extract coefficient and z-term
            match = re.match(r'([+-]?\d*)([*]?)z(\d+)', term)
            if match:
                coeff_str, _, z_num = match.groups()
                # Handle implicit coefficients
                if coeff_str in ('', '+'):
                    coeff = 1
                elif coeff_str == '-':
                    coeff = -1
                else:
                    coeff = int(coeff_str)
                result[f'z{z_num}'] = coeff
            else:
                # Handle cases where the term is just a number (though not expected in this context)
                num_match = re.match(r'([+-]?\d+)', term)
                if num_match:
                    coeff = int(num_match.group(1))
                    # This case is not handled in the result dictionary since there's no z-term
                    pass
        return result


    def latex_cone_tikz(rays):
        angles = get_cummulative_angles(rays)
        index_cones = get_cones(convert_ray_list(rays))
        print(r'\begin{equation}')
        print(r'    \begin{tikzpicture}[scale = 2, baseline=-.5ex]')
        print(r'        \coordinate (A) at (0, 0);')

        for i in range(len(rays)):
            coord = get_coordinate(rays[i])
            print(f'        \coordinate ({alphabet[i + 1]}) at ({coord[0]}, {coord[1]});')

        for i in range(len(rays)):
            print(fr'        \fill [{colors[i]}!60] (A) -- ({alphabet[i + 1]}) arc ({angles[i][0]}:{angles[i][1]}:1.0cm) -- cycle;')

        for i in range(len(rays)):
            print(fr'        \draw[very thick] (A) -- ({alphabet[i + 1]});')
        for i in range(len(rays)):

            mid_point = get_mid_line(index_cones[i])
            print(fr'        \coordinate ({alphabet[len(rays) + 1 + i]}) at ({mid_point[0]}, {mid_point[1]});')
            print(fr'        \draw ({alphabet[len(rays) + 1 + i]}) node[scale=1.5] {{\(\mathbf{{\sigma_{{{i}}}}}\)}};')

        print(r'     \end{tikzpicture}')
        print(r'\end{equation}')
        print("")
        print(r"\newpage")


    def latex_antican(pic):
        string_gen_list = []
        for pic_basis in pic.basis():
            pic_coho = coho(pic_basis.lift())
            temp = str(pic_coho)[1:-1]
            temp = temp.replace(" ", "")
            string_gen_list.append(temp)
        coeffs = []

        for i in range(len(string_gen_list)):
            coeffs.append(extract_coefficients(string_gen_list[i]))

        shape = ((len(fan.rays()) - 2), (len(fan.rays()) - 2))
        antican_matrix = np.zeros(shape)
        for i in range(len(coeffs)):
            for j in range(2, len(fan.rays())):
                try:
                    antican_matrix[j - 2][i] = int(coeffs[i][f"z{j}"])
                except:
                    continue
        antican = 0
        for j in range(len(fan.rays())):
            antican += X.divisor(j)
        antican_coho = str(coho(antican))
        antican_coho = antican_coho.replace("[", "")
        antican_coho = antican_coho.replace("]", "")
        antican_coeffs = extract_coefficients(antican_coho)
        antican_shape = (len(fan.rays()) - 2)
        antican_vector = np.zeros(antican_shape)
        for i in range(2, len(fan.rays())):
            try:
                antican_vector[i - 2] = int(antican_coeffs[f"z{i}"])
            except:
                continue
        antican_solution = np.linalg.solve(antican_matrix, antican_vector)
        antican_solution = np.round(antican_solution).astype(int)
        antican_str = str(coho(antican))[1:-1]
        antican_str = antican_str.replace("*", "")
        antican_str = antican_str.replace(" ", "")
        antican_str = antican_str.replace("z", "z_")

        antican_sol_str = ""
        antican_sol_vector_str = ""
        for i in range(len(antican_solution)):
            if i < len(antican_solution) - 1:
                value = str(antican_solution[i]).replace(".0", "")
                value = value.replace("-0", "0")
                antican_sol_str += f"{value}e_{{{i + 1}}}+"
                antican_sol_vector_str += f"{value},"
            else:
                value = str(antican_solution[i]).replace(".0", "")
                value = value.replace("-0", "0")
                antican_sol_str += f"{value}e_{{{i + 1}}}"
                antican_sol_vector_str += f"{value}"

        str_pic_coho = "    "
        for i in range(len(fan.rays()) - 2):
            if i == 0:
                index = f"e_{{{i + 1}}}&=["
            else:
                index = f"e_{{{i + 1}}}=["
            pic_coho = index + str(coho(pic.gen(i).lift()))[1:-1]
            pic_coho = pic_coho.replace("*", "")
            pic_coho = pic_coho.replace(" ", "")
            pic_coho = re.sub(r'z(\d+)', r'z_{\1}', pic_coho)
            if i < len(fan.rays()) - 3:
                pic_coho += "],"
            else:
                pic_coho += r"]\\"
            str_pic_coho += pic_coho
        print(r"\begin{align*}")
        print(str_pic_coho)
        print(fr"    -\mathcal{{K}}&=\sum_{{i=0}}^{{{len(fan.rays()) - 1}}} D_i=[{antican_str}]\\")
        print(fr"    -\mathcal{{K}}&={antican_sol_str},-\mathcal{{K}}=[{antican_sol_vector_str}]")
        print(r"\end{align*}")
        print("")

        return antican_solution


    def latex_intersection_matrix(intersection_matrix):
        first_div_row = "& "
        for i in range(len(fan.rays())):
            if i < len(fan.rays()) - 1:
                first_div_row += f"D_{{{i}}} & "
            else:
                first_div_row += f"D_{{{i}}}"
        first_div_row += r"\\ \hline"

        c_str = "c" * len(fan.rays())

        print(r"\[")
        print(fr"\begin{{array}}{{|c|{c_str}|}}")
        print(r"\hline")
        print(first_div_row)

        for i, row in enumerate(intersection_matrix):
            row_str = " & ".join(map(str, row))
            print(fr"D_{{{i}}} & {row_str} \\")

        print(r"\hline")
        print(r"\end{array}")
        print(r"\]")
        print("")


    def latex_rescale_antican(antican_solution):
        antican_solution = list(antican_solution)
        max_value = max(antican_solution)
        for i in range(len(antican_solution)):
            if antican_solution[i] == max_value:
                antican_solution[i] = str(1)
                antican_solution[i] += ","
            elif antican_solution[i] == 0:
                antican_solution[i] = f"{0},"
            else:
                antican_solution[i] = fr"\frac{{{str(antican_solution[i])}}}{{{str(max_value)}}},"

        antican_rescale_str = r"Rescaled: $-\mathcal{K}=\left["
        for i in range(len(antican_solution)):
            antican_rescale_str += antican_solution[i]
        antican_rescale_str = antican_rescale_str[:-1]
        antican_rescale_str += r"\right]$"
        print(antican_rescale_str)
        print("")


    def info_rescale_antican(antican_solution):
        antican_solution = list(antican_solution)
        max_coeff = max(antican_solution)
        rescaled_antican = []
        for i in range(len(antican_solution)):
            rescaled = antican_solution[i] / max_coeff
            rescaled_antican.append(rescaled)

        return rescaled_antican


    def latex_Kc_contains(antican_solution):
        print(fr"Kc contains? {Kc.contains(tuple(antican_solution))}")
        print("")
        print("Any area of interest?")
        print("")


    def latex_rays(str_rays):
        print(fr"rays = {str_rays}")
        print("")



    def compute_intersection_matrix(number_of_rays, X):
        '''
        This function computes the intersection matrix of the torus invariant divisors
            corresponding to the rays of the fan, using only Sage

        number_of_rays: is the number of rays of the fan

        X: is the toric variety computed by Sage

        Returns: the intersection matrix of the torus invariant divisors
        '''

        intersection_matrix = [[0 for i in range(number_of_rays)] for j in range(number_of_rays)]
        for i in range(number_of_rays):
            for j in range(number_of_rays):
                i_coho = coho(X.divisor(i))
                j_coho = coho(X.divisor(j))
                intersection_number = X.integrate(i_coho * j_coho)
                intersection_matrix[i][j] = intersection_number

        return intersection_matrix


    alphabet = list(string.ascii_uppercase)
    colors = ['Yellow', 'OrangeRed', 'Rhodamine', 'Violet', 'SkyBlue', 'LimeGreen', 'Blue', 'Red', 'green', 'Orange', 'white']
    for ray_input in input_strings:
        rays_of_fan = []
        for pair in ray_input.split(";"):
            num1, num2 = map(int, [x.strip() for x in pair.split(",")])
            rays_of_fan.append((num1, num2))
        cones_of_fan = define_cones(rays_of_fan)
        fan = create_fan(rays_of_fan, cones_of_fan)

        X = ToricVariety(fan)
        Kc = X.Kaehler_cone()

        pic = X.rational_class_group()
        coho = X.cohomology_ring()
        x_sorted_rays = sort_from_x(ray_input)

        intersection_matrix = compute_intersection_matrix(len(fan.rays()), X)
        file_name = 'CSCKBlHirzebruch25Points'
        print(fr"\textbf{{File name: {file_name}{input_strings.index(ray_input)}}}")
        print("")
        latex_rays(ray_input)
        antican_solution = latex_antican(pic)
        info_antican = []
        for i in range(len(antican_solution)):
            info_antican.append(antican_solution[i])
        info_antican_rescaled = info_rescale_antican(info_antican)
        full_info = np.array([rays_of_fan, info_antican_rescaled], dtype=object)

        np.save(fr'Info/{file_name}Info{input_strings.index(ray_input)}', full_info)

        latex_rescale_antican(antican_solution)
        latex_Kc_contains(antican_solution)
        latex_intersection_matrix(intersection_matrix)
        latex_cone_tikz(x_sorted_rays)

print_everything(input_rays)

# Reset stdout
sys.stdout = sys.__stdout__

# Copy captured output to clipboard
output = output_buffer.getvalue()
pyperclip.copy(output)