import sympy as sp
import numpy as np
from scipy.optimize import linprog


# =========================================================
# SYMBOLS
# =========================================================

def generate_symbolic_coeffs(n, base_name="a"):
    return list(sp.symbols(f"{base_name}0:{n}"))


# =========================================================
# DIVISOR CONSTRUCTION
# =========================================================

def gen_divisor(coeffs, generators):
    return list(np.sum(
        [[c * g for g in generators[i]] for i, c in enumerate(coeffs)],
        axis=0
    ))


# =========================================================
# INTERSECTION THEORY
# =========================================================

def intersection_expressions(surface, coeffs, curves):

    M = sp.Matrix(surface.intersection_matrix)
    D = sp.Matrix(gen_divisor(coeffs, surface.picard_generators))

    results = []

    for c in curves:
        val = (D.T * M * sp.Matrix(c))[0]
        results.append(sp.simplify(val))

    return results


# =========================================================
# LINEAR PROGRAMMING FORM
# =========================================================

def expressions_to_lp(exprs, vars):

    lp = []
    zero = {v: 0 for v in vars}

    for e in exprs:
        e = sp.expand(e)

        coeffs = [e.coeff(v) for v in vars]
        rhs = -e.subs(zero)

        lp.append((coeffs, rhs))

    return lp


def is_redundant(ineqs, i):

    others = [ineqs[j] for j in range(len(ineqs)) if j != i]
    ci, bi = ineqs[i]

    A, b = [], []

    for c, rhs in others:
        A.append([-x for x in c])
        b.append(-rhs)

    A.append(ci)
    b.append(bi - 1e-6)

    res = linprog([0]*len(ci), A_ub=A, b_ub=b, method="highs")

    return not res.success


def remove_redundant(lp):

    return [
        lp[i]
        for i in range(len(lp))
        if not is_redundant(lp, i)
    ]


# =========================================================
# CLEAN FORM CONVERSION
# =========================================================

def move_to_rhs(expr):

    expr = sp.expand(expr)

    pos, neg = [], []

    for t in expr.as_ordered_terms():
        c = t.as_coeff_Mul()[0]

        if c >= 0:
            pos.append(t)
        else:
            neg.append(-t)

    lhs = sp.simplify(sum(pos)) if pos else 0
    rhs = sp.simplify(sum(neg)) if neg else 0

    return lhs, rhs


# =========================================================
# PUBLIC API (USED BY PLOTTER)
# =========================================================

def get_inequalities(surface, curves, base_name="a"):
    """
    Returns:
        list of (lhs, rhs) SymPy expressions
    """

    coeffs = generate_symbolic_coeffs(
        len(surface.picard_generators),
        base_name
    )

    vars = coeffs

    exprs = intersection_expressions(surface, coeffs, curves)

    lp = expressions_to_lp(exprs, vars)
    lp = remove_redundant(lp)

    # keep same indexing but safe filtering
    cleaned = [
        move_to_rhs(exprs[i])
        for i in range(len(exprs))
        if i < len(lp)
    ]

    return cleaned


# =========================================================
# UI FORMATTER (used by Plotter)
# =========================================================

def format_inequalities(surface, curves, latex=True):
    """
    Converts output into display strings.
    """

    ineqs = get_inequalities(surface, curves)

    out = []

    for lhs, rhs in ineqs:

        if latex:
            out.append(f"${sp.latex(lhs)} > {sp.latex(rhs)}$")
        else:
            out.append(f"{lhs} > {rhs}")

    return out