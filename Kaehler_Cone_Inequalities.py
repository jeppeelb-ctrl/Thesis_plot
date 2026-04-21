import numpy as np
import itertools
import sympy as sp
import numpy as np
from scipy.optimize import linprog


def expressions_to_lp(expressions, variables):
    """
    Convert list of SymPy expressions (all assumed > 0)
    into LP format: (coeff_vector, rhs)
    """

    lp_form = []

    for expr in expressions:
        expr = sp.expand(expr)

        coeffs = [sp.expand(expr).coeff(v) for v in variables]

        # constant term
        rhs = -expr.subs({v: 0 for v in variables})

        lp_form.append((coeffs, rhs))

    return lp_form


def is_redundant(ineqs, i):
    """
    Check if inequality i is redundant.
    ineqs: list of (coeff_vector, rhs) representing coeff·x >= rhs
    """

    # Remove i-th inequality
    others = [ineqs[j] for j in range(len(ineqs)) if j != i]

    # We try to violate inequality i:
    # coeff_i · x < rhs_i  →  coeff_i · x <= rhs_i - ε
    coeff_i, rhs_i = ineqs[i]

    # Convert to standard form: A x <= b
    A = []
    b = []

    for coeff, rhs in others:
        A.append([-c for c in coeff])   # ≥ → ≤
        b.append(-rhs)

    # Add violation condition
    epsilon = 1e-6
    A.append(coeff_i)
    b.append(rhs_i - epsilon)

    # Objective doesn't matter
    c = [0]*len(coeff_i)

    res = linprog(c, A_ub=A, b_ub=b, method='highs')

    # If infeasible → cannot violate → redundant
    return not res.success


def find_redundancies(ineqs):
    redundant = []
    for i in range(len(ineqs)):
        if is_redundant(ineqs, i):
            redundant.append(i)
    return redundant


def generate_symbolic_coeffs(n, base_name="a"):
    """
    ստեղծ [a0, a1, ..., a_{n-1}] or [a, b, c, ...] if n <= 26
    """
    if base_name == "letters" and n <= 26:
        letters = "abcdefghijklmnopqrstuvwxyz"
        return list(sp.symbols(" ".join(letters[:n])))
    else:
        return list(sp.symbols(f"{base_name}0:{n}"))


def symbolic_picard_divisor(Pic_generators, base_name="a"):
    """
    Returns ONE symbolic divisor in Picard basis:
    replaces coefficients by variables and applies gen_Pic_divisor
    """
    coeffs = generate_symbolic_coeffs(len(Pic_generators), base_name)
    return coeffs, gen_Pic_divisor(coeffs, Pic_generators)


def all_symbolic_nef_divisors(Pic_generators, base_name="a"):
    """
    Returns symbolic expressions for ALL Kaehler nef check curves,
    but expressed in Picard basis.
    """
    # original binary curves
    curves = Kaehler_check_curves(len(Pic_generators))

    # symbolic variables
    coeffs = generate_symbolic_coeffs(len(Pic_generators), base_name)

    symbolic_results = []

    for curve in curves:
        # replace 1 -> variable, 0 -> 0
        symbolic_curve = [
            coeffs[i] if curve[i] == 1 else 0
            for i in range(len(curve))
        ]

        divisor = gen_Pic_divisor(symbolic_curve, Pic_generators)
        symbolic_results.append(divisor)

    return coeffs, symbolic_results

def intersection_calculator(v1_vec, v2_vec, inter_matrix):

    v2_transpose = np.transpose(v2_vec)
    intermediate_result = np.dot(v1_vec, inter_matrix)

    return np.dot(intermediate_result, v2_transpose)

def Kaehler_check_curves(Picard_rank):
    return [list(t) for t in itertools.product([0, 1], repeat=Picard_rank)]


def print_intersection_matrix(inter_matrix):

    temp_str = "    "
    for j in range(len(inter_matrix)):
        temp_str += "D" + str(j) + " "
    print(temp_str)
    for j in range(len(inter_matrix)):
        print("D" + str(j), inter_matrix[j])


def sum_coeff_list(lists):
    return list(np.sum(lists, axis=0))


def compute_intersection_matrix(number_of_rays, rays):

    intersection_matrix = [[0 for i in range(number_of_rays)] for j in range(number_of_rays)]
    for i in range(number_of_rays):
        if i == number_of_rays - 1:
            left_ray = rays[i - 1]
            right_ray = rays[0]
        else:
            left_ray = list(rays[i - 1])
            right_ray = list(rays[i + 1])
        intersection_number = 0
        for number in range(-15, 15):
            added_rays = [left_ray[0] + right_ray[0], left_ray[1] + right_ray[1]]
            numb_rays = [number * rays[i][0], number * rays[i][1]]
            if added_rays == numb_rays:
                intersection_number = -number
        intersection_matrix[i][i] = intersection_number
        if i < number_of_rays - 1:
            intersection_matrix[i][i + 1] = 1
            intersection_matrix[i + 1][i] = 1

    intersection_matrix[0][number_of_rays - 1] = 1
    intersection_matrix[number_of_rays - 1][0] = 1

    return intersection_matrix


def Picard_generators(rays_of_fan):
    D_description = []
    D_zero_description = [0, 0]
    D_one_description = [0, 0]

    for i in range(2, len(rays_of_fan)):
        one_description = - int(np.inner([1, 0], rays_of_fan[i]))
        D_one_description.append(one_description)

        zero_description = int(np.inner([rays_of_fan[1][1], -1], rays_of_fan[i]))
        D_zero_description.append(zero_description)

    D_description.append(D_zero_description)
    D_description.append(D_one_description)

    for j in range(2, len(rays_of_fan) - 2):
        description = [0 for i in range(len(rays_of_fan))]
        description[j] = 1
        D_description.append(description)

    return D_description


def gen_Pic_divisor(coeffs, Pic_generators):
    tmp_div = []

    for i in range(len(coeffs)):
        tmp = []
        for k in range(len(Pic_generators[i])):
            div_Pic_gen_coeff = coeffs[i] * Pic_generators[i][k]
            tmp.append(div_Pic_gen_coeff)
        tmp_div.append(tmp)

    div = sum_coeff_list(tmp_div)
#    for _ in range(len(div)):
#        div[_] = int(div[_])
    return div


def intersection_all_nef(symbolic_coeffs, Pic_generators, nef_curves, intersection_matrix):
    """
    Computes intersection of ONE symbolic divisor with ALL nef curves
    """
    im = sp.Matrix(intersection_matrix)

    # symbolic divisor in toric divisor basis
    sym_div = sp.Matrix(gen_Pic_divisor(symbolic_coeffs, Pic_generators))

    results = []

    for curve in nef_curves:
        curve_vec = sp.Matrix(curve)
        val = (sym_div.T * im * curve_vec)[0]
        results.append(sp.simplify(val))

    return sym_div, results


info = np.load(f"CasesOfInterestInfo\CSCKBlP1xP11PointsInfo1.npy", allow_pickle=True)
rays = info[0]

Pic_gen = Picard_generators(rays)

print(f"Initial rays = {rays}")
Kaehler_nef_check_curves = Kaehler_check_curves(len(Pic_gen))

nef_check_curves = []
for i in range(len(Kaehler_nef_check_curves)):
    nef_check_curves.append(gen_Pic_divisor(Kaehler_nef_check_curves[i], Pic_gen))

intersection_matrix = compute_intersection_matrix(len(rays), rays)
im = np.array(intersection_matrix)
#print_intersection_matrix(intersection_matrix)

coeffs, div = symbolic_picard_divisor(Pic_gen, base_name="letters")

sym_div, intersections = intersection_all_nef(
    coeffs,
    Pic_gen,
    nef_check_curves,
    intersection_matrix
)

intersections = intersections[1:]

vars = generate_symbolic_coeffs(len(Pic_gen), base_name="letters")

lp = expressions_to_lp(intersections, vars)

redundant = set(find_redundancies(lp))

non_redundant = [
    intersections[i]
    for i in range(len(intersections))
    if i not in redundant
]

for ineq in non_redundant:
    print(str(ineq) + ' > 0')
