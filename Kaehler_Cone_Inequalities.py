import sympy as sp
import numpy as np
from scipy.optimize import linprog
from dataclasses import dataclass


# =========================
# Data structures
# =========================

@dataclass
class ToricSurface:
    rays: np.ndarray
    picard_generators: list
    intersection_matrix: np.ndarray


# =========================
# Utility functions
# =========================

def generate_symbolic_coeffs(n, base_name="a"):
    if base_name == "letters" and n <= 26:
        letters = "abcdefghijklmnopqrstuvwxyz"
        return list(sp.symbols(" ".join(letters[:n])))
    return list(sp.symbols(f"{base_name}0:{n}"))


def sum_coeff_list(lists):
    return list(np.sum(lists, axis=0))


# =========================
# Geometry construction
# =========================

def compute_intersection_matrix(number_of_rays, rays):
    SEARCH_RANGE = range(-15, 15)

    intersection_matrix = [[0]*number_of_rays for _ in range(number_of_rays)]

    for i in range(number_of_rays):
        left_ray = rays[i - 1]
        right_ray = rays[(i + 1) % number_of_rays]

        intersection_number = 0
        for number in SEARCH_RANGE:
            if [
                left_ray[0] + right_ray[0],
                left_ray[1] + right_ray[1]
            ] == [number * rays[i][0], number * rays[i][1]]:
                intersection_number = -number

        intersection_matrix[i][i] = intersection_number

        intersection_matrix[i][(i + 1) % number_of_rays] = 1
        intersection_matrix[(i + 1) % number_of_rays][i] = 1

    return np.array(intersection_matrix)


def Picard_generators(rays):
    D_description = []

    D_zero = [0, 0]
    D_one = [0, 0]

    for i in range(2, len(rays)):
        D_one.append(-int(np.inner([1, 0], rays[i])))
        D_zero.append(int(np.inner([rays[1][1], -1], rays[i])))

    D_description.append(D_zero)
    D_description.append(D_one)

    for j in range(2, len(rays) - 2):
        desc = [0]*len(rays)
        desc[j] = 1
        D_description.append(desc)

    return D_description


def build_surface(rays):
    pic = Picard_generators(rays)
    inter = compute_intersection_matrix(len(rays), rays)
    return ToricSurface(rays=rays, picard_generators=pic, intersection_matrix=inter)


# =========================
# Divisors & intersections
# =========================

def gen_Pic_divisor(coeffs, Pic_generators):
    tmp = []
    for i in range(len(coeffs)):
        tmp.append([coeffs[i] * g for g in Pic_generators[i]])
    return sum_coeff_list(tmp)


def symbolic_picard_divisor(Pic_generators, base_name="a"):
    coeffs = generate_symbolic_coeffs(len(Pic_generators), base_name)
    return coeffs, gen_Pic_divisor(coeffs, Pic_generators)


def tmp_Kaehler_test_curves(rank):
    tmp = [[0]*(rank+2) for _ in range(rank+2)]
    for i in range(rank+2):
        tmp[i][i] = 1
    return tmp


def intersection_all_nef(surface, coeffs, nef_curves):
    im = sp.Matrix(surface.intersection_matrix)
    sym_div = sp.Matrix(gen_Pic_divisor(coeffs, surface.picard_generators))

    results = []
    for curve in nef_curves:
        val = (sym_div.T * im * sp.Matrix(curve))[0]
        results.append(sp.simplify(val))

    return sym_div, results


# =========================
# LP + redundancy
# =========================

def expressions_to_lp(expressions, variables):
    lp_form = []

    zero_subs = {v: 0 for v in variables}

    for expr in expressions:
        expr = sp.expand(expr)
        coeffs = [expr.coeff(v) for v in variables]
        rhs = -expr.subs(zero_subs)

        lp_form.append((coeffs, rhs))

    return lp_form


def is_redundant(ineqs, i):
    others = [ineqs[j] for j in range(len(ineqs)) if j != i]

    coeff_i, rhs_i = ineqs[i]

    A = []
    b = []

    for coeff, rhs in others:
        A.append([-c for c in coeff])
        b.append(-rhs)

    epsilon = 1e-6
    A.append(coeff_i)
    b.append(rhs_i - epsilon)

    c = [0]*len(coeff_i)

    res = linprog(c, A_ub=A, b_ub=b, method='highs')

    return not res.success


def find_redundancies(ineqs):
    return [i for i in range(len(ineqs)) if is_redundant(ineqs, i)]


def move_negatives_to_rhs(expr):
    """
    Takes a SymPy expression f (interpreted as f > 0)
    and returns (lhs, rhs) such that:
        lhs > rhs
    with no negative coefficients.
    """

    expr = sp.expand(expr)

    pos_terms = []
    neg_terms = []

    for term in expr.as_ordered_terms():
        coeff = term.as_coeff_Mul()[0]

        if coeff >= 0:
            pos_terms.append(term)
        else:
            neg_terms.append(-term)  # move to RHS, flip sign

    lhs = sp.simplify(sum(pos_terms)) if pos_terms else sp.Integer(0)
    rhs = sp.simplify(sum(neg_terms)) if neg_terms else sp.Integer(0)

    return lhs, rhs


# =========================
# Pipeline
# =========================

def compute_nonredundant_inequalities(surface):
    coeffs, _ = symbolic_picard_divisor(surface.picard_generators, base_name="letters")

    nef_curves = tmp_Kaehler_test_curves(len(surface.picard_generators))

    _, intersections = intersection_all_nef(surface, coeffs, nef_curves)

    vars = generate_symbolic_coeffs(len(surface.picard_generators), base_name="letters")

    lp = expressions_to_lp(intersections, vars)
    redundant = set(find_redundancies(lp))

    non_redundant =  [
        intersections[i]
        for i in range(len(intersections))
        if i not in redundant
    ]

    return [
        move_negatives_to_rhs(ineq)
        for ineq in non_redundant
    ]


def get_inequalities(inequalities):
    ineqs = []
    for (lhs, rhs) in inequalities:
        inequality = f"{lhs} > {rhs}"
        ineqs.append(inequality)
    return ineqs

# =========================
# Main
# =========================

def main():
    info = np.load("CasesOfInterestInfo/CSCKBlP22PointsInfo8.npy", allow_pickle=True)
    rays = info[0]

    surface = build_surface(rays)

    print(f"Initial rays = {rays}")

    inequalities = compute_nonredundant_inequalities(surface)

    #for (lhs, rhs) in inequalities:
    #    print(f"{lhs} > {rhs}")
