import string
import math
from math import sqrt
import numpy as np
import sys
import pyperclip
from io import StringIO

color_cases = {
    8: [
        'RedOrange!60', 'Dandelion!60', 'SpringGreen!60',
        'YellowGreen!60', 'SkyBlue!60', 'CornflowerBlue!60',
        'Mulberry!60', 'Orchid!60'
    ],
    7: [
        'RedOrange!60', 'Dandelion!60', 'SpringGreen!60',
        'YellowGreen!60', 'SkyBlue!60', 'CornflowerBlue!60',
        'Mulberry!60'
    ],
    6: [
        'RedOrange!60', 'Dandelion!60', 'SpringGreen!60',
        'YellowGreen!60', 'CornflowerBlue!60', 'Mulberry!60'
    ],
    5: [
        'RedOrange!60', 'Dandelion!60',
        'YellowGreen!60', 'CornflowerBlue!60', 'Mulberry!60'
    ],
    4: [
        'RedOrange!60', 'YellowGreen!60',
        'CornflowerBlue!60', 'Mulberry!60'
    ],
    3: [
        'RedOrange!60', 'YellowGreen!60', 'CornflowerBlue!60'
    ]
}

# Redirect stdout to capture print statements
output_buffer = StringIO()
sys.stdout = output_buffer

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

    fan_str_rays = []
    for ray in sorted_points:
        ray = list(ray)
        fan_str_rays.append(ray)
    return fan_str_rays


def convert_to_jupyter_ray_input(rays):
    str_rays = str(rays)
    str_rays = str_rays.replace("[","")
    str_rays = str_rays.replace("]", "")
    str_rays = str_rays.replace("(", "")
    str_rays = str_rays.replace(")", "")
    str_rays = str_rays.replace(" ", "")

    elements = str_rays.split(",")
    result = elements[0]
    for i in range(1, len(elements)):
        if i % 2 == 0:
            result += ";" + elements[i]
        else:
            result += "," + elements[i]

    return result


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


def latex_cone_tikz(rays):
    angles = get_cummulative_angles(rays)
    index_cones = get_cones(convert_ray_list(rays))
    print(r'\begin{equation*}')
    print(r'    \begin{tikzpicture}[scale = 2, baseline=-.5ex]')
    print(r'        \coordinate (A) at (0, 0);')

    for i in range(len(rays)):
        coord = get_coordinate(rays[i])
        print(f'        \coordinate ({alphabet[i + 1]}) at ({coord[0]}, {coord[1]});')

    for i in range(len(rays)):
        print(fr'        \fill [{color_cases.get(len(rays))[i]}] (A) -- ({alphabet[i + 1]}) arc ({angles[i][0]}:{angles[i][1]}:1.0cm) -- cycle;')

    for i in range(len(rays)):
        print(fr'        \draw[very thick] (A) -- ({alphabet[i + 1]});')
    for i in range(len(rays)):

        mid_point = get_mid_line(index_cones[i])
        print(fr'        \coordinate ({alphabet[len(rays) + 1 + i]}) at ({mid_point[0]}, {mid_point[1]});')
        print(fr'        \draw ({alphabet[len(rays) + 1 + i]}) node[scale=1.5] {{\(\mathbf{{\sigma_{{{i}}}}}\)}};')

    print(r'     \end{tikzpicture}')
    print(r'\end{equation*}')
    print("")


alphabet = list(string.ascii_uppercase)


info = np.load(f"CasesOfInterestInfo\CSCKBlP22PointsInfo1.npy", allow_pickle=True)

initial_rays = info[0]

sorted_rays = sort_from_x(convert_to_jupyter_ray_input(initial_rays))

latex_cone_tikz(sorted_rays)

# Reset stdout
sys.stdout = sys.__stdout__

# Copy captured output to clipboard
output = output_buffer.getvalue()
pyperclip.copy(output)
