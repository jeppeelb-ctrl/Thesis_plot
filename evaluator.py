# evaluator.py
import numpy as np

# =========================
# Utility functions
# =========================

def sum_coeff_list(lists):
    return list(np.sum(lists, axis=0))


def check_non_negative(lst):

    if any(x <= 0 for x in lst):
        return False
    return True


def split_theta(theta):
    split_index = int(len(theta) / 2)
    alpha = theta[:split_index]
    beta = theta[split_index:]
    return alpha, beta


# =========================
# Divisors & intersections
# =========================

def gen_Pic_divisor(coeffs, Pic_generators):
    tmp = []
    for i in range(len(coeffs)):
        tmp.append([coeffs[i] * g for g in Pic_generators[i]])
    return sum_coeff_list(tmp)


def generate_invariant_divisor_i(number_of_rays, index):
    invariant_div = [0 for _ in range(number_of_rays)]
    for i in range(number_of_rays):
        if i == index:
            invariant_div[i] = 1

    return invariant_div


def convert_to_pic_description(div, surface):
    pic_gen = surface.picard_generators
    converted_div = [0 for _ in range(len(pic_gen[0]))]
    for i, coeff in enumerate(div):
        for j, val in enumerate(pic_gen[i]):
            converted_div[j] += coeff * val
    return converted_div


# =========================
# Intersection number calculations
# =========================

def intersection_calculator(div1, div2, surface):
    im = surface.intersection_matrix
    if len(div1) == len(surface.picard_generators):
        conv_div1 = convert_to_pic_description(div1, surface)
    else:
        conv_div1 = np.array(div1)
    if len(div2) == len(surface.picard_generators):
        conv_div2 = convert_to_pic_description(div2, surface)
    else:
        conv_div2 = np.array(div2)

    div2_transpose = np.transpose(conv_div2)
    intermediate_result = np.dot(conv_div1, im)

    return np.dot(intermediate_result, div2_transpose)


def is_kaehler(divisor, surface):
    intersection_numbers = []
    for i in range(len(surface.rays)):
        invariant_divisor = generate_invariant_divisor_i(len(surface.rays), i)
        intersection_number = intersection_calculator(div1=divisor, div2=invariant_divisor, surface=surface)
        intersection_numbers.append(intersection_number)

    return check_non_negative(intersection_numbers), intersection_numbers


def compute_anti_canonical(surface):
    pic_gen = surface.picard_generators
    anti_canonical = [1 for _ in range(len(pic_gen) + 2)]
    #anti_canonical = convert_to_pic_description(sum_of_divs, surface)
    return anti_canonical

def compute_alpha_squared(alpha, surface):
    alpha_squared = intersection_calculator(div1=alpha, div2=alpha, surface=surface)
    return alpha_squared


def compute_alpha_beta(alpha, beta, surface):
    alpha_beta = intersection_calculator(div1=alpha, div2=beta, surface=surface)
    return alpha_beta


def compute_full_factor_Jeq(alpha, beta, surface):

    result = np.inf
    invariant_divisor = 0
    is_alpha_kaehler, alpha_invariant_intersections = is_kaehler(alpha, surface)
    is_beta_kaehler, beta_invariant_intersections = is_kaehler(beta, surface)
    if is_alpha_kaehler and is_beta_kaehler:
        alpha_squared = compute_alpha_squared(alpha=alpha, surface=surface)
        alpha_beta = compute_alpha_beta(alpha=alpha, beta=beta, surface=surface)
        first_factor = 2 * (alpha_beta / alpha_squared)
        for i in range(len(alpha_invariant_intersections)):
            second_factor = (beta_invariant_intersections[i] / alpha_invariant_intersections[i])
            full_factor = first_factor - second_factor
            if full_factor < result:
                result = full_factor
                invariant_divisor = i
    else:

        return np.nan, False
    return result, True


def compute_full_factor_cscK(alpha, beta, surface):
    result = -np.inf
    antican = compute_anti_canonical(surface)
    invariant_divisor = 0
    is_alpha_kaehler, alpha_invariant_intersections = is_kaehler(alpha, surface)
    _, antican_invariant_intersections = is_kaehler(antican, surface)
    if is_alpha_kaehler:
        alpha_squared = compute_alpha_squared(alpha=alpha, surface=surface)
        alpha_beta = compute_alpha_beta(alpha=alpha, beta=antican, surface=surface)
        first_factor = 2 * (alpha_beta / alpha_squared)
        for i in range(len(alpha_invariant_intersections)):
            second_factor = (antican_invariant_intersections[i] / alpha_invariant_intersections[i])
            full_factor = first_factor - second_factor
            if full_factor > result:
                result = full_factor
                invariant_divisor = i
    else:

        return np.nan, False
    return -result, True

EQUATIONS = {
    "J(α,β)": compute_full_factor_Jeq,
    "I(α)": compute_full_factor_cscK,
}


def evaluate_point(name, theta, surface):
    fn = EQUATIONS.get(name)
    if fn is None:
        raise ValueError(f"Unknown equation: {name!r}. Available: {list(EQUATIONS)}")
    theta = np.asarray(theta)
    divisors = split_theta(theta)
    alpha = divisors[0][0]
    beta = divisors[1][0]
    return fn(alpha, beta, surface)
