import argparse
import os
import re
import math
import numpy as np
from sage.all import *


def intersection_calculator(v1_coeff, v2_coeff, intersection_matrix, fan, gen_list_count, gen_list):
    '''
    This function calculates the intersection of two divisors using vectors and matrices

    v1_coeff: is the coefficients of the first divisor
    v2_coeff: is the coefficients of the second divisor

    Returns the intersection number of v1_coeff and v2_coeff
    '''

    v1_vec = convert_to_coho(v1_coeff, gen_list_count, fan)
    if len(v2_coeff) == len(gen_list):
        v2_vec = np.array(convert_to_coho(v2_coeff, gen_list_count, fan))
    else:
        v2_vec = np.array(v2_coeff)

    v2_transpose = np.transpose(v2_vec)
    intermediate_result = np.dot(v1_vec, intersection_matrix)

    return np.dot(intermediate_result, v2_transpose)


def extract_coefficients(s):
    '''
    This function uses regular expressions to extract the information computed by Sagemath
        of the Picard group.

    s: is a string of containing the information computed by Sagemath of a generator of the Picard group.

    Returns: 'result' a list containing the coefficients values, represented by the torus invariant divisors.
    '''
    # Split the string into terms
    terms = re.findall(r'([+-]?\d*[*]?z\d+|[+-]?z\d+)', s)

    result = {}
    for term in terms:
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
    return result


def generate_pic_gen_lists(picard_group, fan, coho):
    '''
    This function creates the generator list of the Picard group, which is the divisor classes [1,0,...]
        [0,1,...] and so forth. Furthermore, it computes a list of the generators in its
        representation of the invariant torus divisors corresponding to the rays of the fan.

    picard_group: is the Picard group computed by Sage

    Returns: 'gen_list' which is a list containing the generators of the Picard group
             'gen_list_count' which is a list containing the information about
                 the generators of the Picard group in its representation using
                 the torus invariant divisors corresponding to the rays
    '''

    gen_list = []
    string_gen_list = []
    for pic_basis in picard_group.basis():
        gen_list.append(pic_basis)
        pic_coho = coho(pic_basis.lift())
        temp = str(pic_coho)[1:-1]
        temp = temp.replace(" ", "")
        string_gen_list.append(temp)

    gen_list_count = [[0 for _ in range(2 * len(fan.rays()))] for _ in range(len(fan.rays()) - 2)]

    coeffs = []

    for i in range(len(string_gen_list)):
        coeffs.append(extract_coefficients(string_gen_list[i]))

    for j in range(len(coeffs)):
        for k in range(len(fan.rays())):
            try:
                if coeffs[j][f"z{k}"] < 0:
                    gen_list_count[j][k + len(fan.rays())] = coeffs[j][f"z{k}"]
                else:
                    gen_list_count[j][k] = coeffs[j][f"z{k}"]
            except:
                continue
    return gen_list, gen_list_count



def first_factor(alpha, beta, intersection_matrix, coho, gen_list_count, gen_list):
    '''
    This function computes the factor 2·(α.β)/α^2

    alpha: is the divisor corresponding to the alpha part of the factor

    beta: is the divisor corresponding to the beta part of the factor

    Returns: the result of the above mentioned computation
    '''

    numerator = intersection_calculator(tuple(alpha), tuple(beta), intersection_matrix, coho, gen_list_count, gen_list)
    denominator = intersection_calculator(tuple(alpha), tuple(alpha), intersection_matrix, coho, gen_list_count, gen_list)
    if denominator == 0:
        return 0
    else:
        result = 2 * (numerator / denominator)

    return result


def second_factor(alpha, beta, C, intersection_matrix, coho, gen_list_count, gen_list):
    '''
    This function computes the factor (β.C)/(α.C)

    alpha: is the divisor corresponding to the alpha part of the factor

    beta: is the divisor corresponding to the beta part of the factor

    Returns: the result of the computation
    '''

    numerator = intersection_calculator(tuple(beta), tuple(C), intersection_matrix, coho, gen_list_count, gen_list)
    denominator = intersection_calculator(tuple(alpha), tuple(C), intersection_matrix, coho, gen_list_count, gen_list)
    if denominator == 0:
        return None
    else:
        result = numerator / denominator

    return result


def convert_to_coho(coefficients, gen_list_count, fan):
    '''
    This function transforms a divisor from the representation in the Picard group
        to the representation using the torus invariant divisors corresponding
        to the rays of the fan

    coefficients: is a tuple of the coefficients of the representation in the Picard group

    gen_list_count: is a list of the generators of the Picard group in the
                    representation using the torus invariant divisors
                    corresponding to the rays of the fan
    '''

    coho_divisor = np.zeros(len(fan.rays()))
    for j in range(len(gen_list_count)):
        for i in range(len(fan.rays())):
            coho_divisor[i] += gen_list_count[j][i] * coefficients[j]
            coho_divisor[i] += gen_list_count[j][i + len(fan.rays())] * coefficients[j]

    return coho_divisor


def convert_C_to_ray_divisor(C):
    '''
    This function converts a torus invariant divisor from its representation in Sage
        to a vector of the form [1,0,...] with 1 in the index of the corresponding ray

    C: is a torus invariant divisor
    '''

    C_div = ''
    for i in range(len(C)):
        if C[i] == 1:
            C_div = ''.join(str(i))
        else:
            continue

    return C_div


def get_coho_C_divisors(number_of_rays):
    '''
    This function creates a list of the torus invariant divisors of the form [1,0...] etc.

    number_of_rays: is the number of rays of the fan

    Returns: a list with every torus invariant divisor of the form [1,0...] etc.
    '''

    C_div_list = [[0 for _ in range(number_of_rays)] for _ in range(number_of_rays)]
    for i in range(len(C_div_list)):
        C_div_list[i][i] = 1

    return C_div_list


def fixed_value_to_index(value, size_of_matrix):
    '''
    This function transforms the fixed value of a coefficient to the closest index in the matrix
    depending on the size of the matrix.
    For example a fixed value of '0.5' for a matrix of size 11 gives the index '5'.

    value: is the value of the fixed value.

    Returns the index in the matrix corresponding to the fixed value.
    '''

    temp = np.linspace(0, 1, size_of_matrix)
    temp2 = np.zeros(len(temp))
    for i in range(len(temp)):
        dist = abs(value - temp[i])
        temp2[i] = dist

    return np.argmin(temp2)


def optim_heat_map_Jeq(alpha, beta, intersection_matrix, C_div_list, fan, gen_list_count, gen_list):
    '''
    This function computes the minimum value of the J-eq quantity J(α,β) checking for
        each torus invariant divisor

    alpha: is the divisor corresponding to the alpha part of the quantity J(α,β)

    beta: is the divisor corresponding to the beta part of the quantity J(α,β)

    Returns: a list with the minimum value at index 0,
             and the torus invariant divisor that realised the minimum value at index 1
    '''

    div_by_zero = []
    value = float(50000)
    C_div = []
    alpha_beta_part = first_factor(alpha, beta, intersection_matrix, fan, gen_list_count, gen_list)
    for C in C_div_list:
        C_part = second_factor(alpha, beta, C, intersection_matrix, fan, gen_list_count, gen_list)
        if C_part is None:
            div_by_zero.append(C)
        else:
            C_value = alpha_beta_part - C_part
            if C_value < value:
                value = C_value
                C_div = C

    return [value, C_div]


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


def compute_intersection_matrix(number_of_rays, X, coho):
    '''
    This function computes the intersection matrix of the torus invariant divisors
        corresponding to the rays of the fan, using only Sage

    number_of_rays: is the number of rays of the fan

    X: is the toric variety computed by Sage

    Returns: the intersection matrix of the torus invariant divisors
    '''

    intersection_matrix = [[0 for _ in range(number_of_rays)] for _ in range(number_of_rays)]
    for i in range(number_of_rays):
        for j in range(number_of_rays):
            i_coho = coho(X.divisor(i))
            j_coho = coho(X.divisor(j))
            intersection_number = X.integrate(i_coho * j_coho)
            intersection_matrix[i][j] = intersection_number

    return intersection_matrix


def process_alpha(alpha_indices, beta, intersection_matrix, Kc, C_div_list, fan, gen_list_count, gen_list, size_of_matrix):
    '''
    This function initiates the computational functions

    Params: the needed information in order to compute the quantity I(α)

    Returns: a tuple with the coefficient values of α, the value of the quantity I(α), and the torus invariant
                divisor that realized the maximum
    '''
    alpha = np.array(alpha_indices)
    if Kc.contains(tuple(alpha)):
        result = optim_heat_map_Jeq(alpha, beta, intersection_matrix, C_div_list, fan, gen_list_count, gen_list)
        value = result[0]
        C_div = convert_C_to_ray_divisor(result[1])
        indices = [fixed_value_to_index(alpha[k], size_of_matrix) for k in range(len(alpha))]
        return (tuple(indices), value, C_div)
    else:
        indices = [fixed_value_to_index(alpha[k], size_of_matrix) for k in range(len(alpha))]
        return (tuple(indices), 50000, None)


def numpy_cartesian_product(arrays):
    '''
    This function produces a meshgrid of the coefficient values of α

    arrays: a list of linspaces, representing the coefficient values of α

    Returns: a column stack with the coefficient values of α
    '''
    mesh = np.meshgrid(*arrays, indexing='ij')
    return np.column_stack([m.ravel() for m in mesh])


def compute_data(beta, basename, output_dir, intersection_matrix, fan, Kc, C_div_list, coho, gen_list_count, gen_list, size_of_matrix):
    '''
    This function initiates the matrices for the computed data, and calls the computational functions
        lastly it saves the computed data

    Params: the needed information in order to determine the size of the matrices, and to compute the quantity I(α)

    '''
    alpha_coeffs = [np.linspace(0, 1, size_of_matrix) for _ in range(len(fan.rays()) - 2)]

    shape = (size_of_matrix,) * (len(fan.rays()) - 2)
    sign_matrix = np.zeros(shape)
    div_matrix = np.empty(shape, dtype=object)

    for alpha_indices in numpy_cartesian_product(alpha_coeffs):
        indices, value, C_div = process_alpha(tuple(alpha_indices), beta, intersection_matrix, Kc, C_div_list, fan, gen_list_count, gen_list, size_of_matrix)

        sign_matrix[indices] = value
        div_matrix[indices] = C_div

    np.save(os.path.join(output_dir, f'{basename}.SignData.npy'), sign_matrix)
    np.save(os.path.join(output_dir, f'{basename}.DivData.npy'), div_matrix)


def find_matrix_size(dim):
    '''
    This function determines the size of the matrices, such that we compute roughly 10^7 points

    dim: the dimension of the Picard group

    Returns: the size of the matrices, in order to compute roughly 10^7 points
    '''
    matrix_size = pow(10, 7 / dim)

    return math.floor(matrix_size)


def compute_with_data(data, output_dir, basename):
    '''
    This function initiates the computation, if the info file is of the correct type

    data: the initial data info, rays of the fan and the values of the fixed divisor
    output_dir: the directory where the info is stored, such that we save the data the same place
    basename: the name of the info file, such that we can efficiently differentiate the computed data

    '''
    # Example computation: replace with your actual computations
    print(f"Performing computations with {basename} \n") # dlf
    if isinstance(data, np.ndarray):

        rays = data[0]
        beta = data[1]
        cones = define_cones(rays)
        fan = Fan(cones=cones, rays=rays)
        number_of_rays = len(fan.rays())
        X = ToricVariety(fan)
        Kc = X.Kaehler_cone()
        pic = X.rational_class_group()
        coho = X.cohomology_ring()
        size_of_matrix = find_matrix_size(number_of_rays - 2)


        gen_list, gen_list_count = generate_pic_gen_lists(pic, fan, coho)
        C_div_list = get_coho_C_divisors(len(fan.rays()))

        intersection_matrix = compute_intersection_matrix(len(fan.rays()), X, coho)

        compute_data(beta, basename, output_dir, intersection_matrix, fan, Kc, C_div_list, coho, gen_list_count, gen_list, size_of_matrix)
    else:
        print("\n\nData type:", type(data))
        print("\n\nData content:", data)


def main():
    parser = argparse.ArgumentParser(description="Load a file and extract folder and basename.")
    parser.add_argument("-file", required=True, help="Path to the file")

    args = parser.parse_args()

    # Extract folder and basename
    file_path = os.path.abspath(args.file)
    folder = os.path.dirname(file_path)
    basename = os.path.basename(file_path)

    #print(f"Full path: {file_path}") # dlf
    #print(f"Folder: {folder}") #dlf
    #print(f"Basename: {basename}") #dlf

    if args.file:
        if os.path.isfile(args.file):
            print(f"Loading file: {args.file}")

            # Extract the directory from the file path
            output_dir = os.path.dirname(args.file)

            # Load the file based on its type
            if args.file.endswith('.npy'):
                data = np.load(args.file, allow_pickle=True)
            elif args.file.endswith('.txt'):
                with open(args.file, 'r') as f:
                    data = f.read()
            else:
                print(f"Unsupported file type: {args.file}")
                return
            print(f"\n{basename} loaded successfully.") #dlf

            # Call the function to compute with the loaded data and save results to the same directory
            compute_with_data(data, output_dir, basename)
        else:
            print(f"File does not exist: {args.file}")
    else:
        print("No file specified.")

if __name__ == "__main__":
    main()
