# geometry.py
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from dataclasses import dataclass


@dataclass
class ToricSurface:
    rays: np.ndarray
    picard_generators: list
    intersection_matrix: np.ndarray
    divisor_dimension: int


class GeometryBuilder:
    """
    Builds toric surface from rays.
    """

    def __init__(self, rays):
        self.rays = np.array(rays)
        self.n = len(rays)

    def build(self):
        pic = self.compute_picard_generators()
        inter = self.compute_intersection_matrix()
        dimension = len(pic)
        return ToricSurface(
            rays=self.rays,
            picard_generators=pic,
            intersection_matrix=inter,
            divisor_dimension=dimension
        )


    # -------------------------
    # Geometry logic
    # -------------------------

    def compute_intersection_matrix(self):
        SEARCH_RANGE = range(-15, 15)

        n = self.n
        rays = self.rays

        M = np.zeros((n, n), dtype=int)

        for i in range(n):
            left = rays[i - 1]
            right = rays[(i + 1) % n]

            intersection_number = 0

            for k in SEARCH_RANGE:
                if np.array_equal(
                        left + right,
                        k * rays[i]
                ):
                    intersection_number = -k

            M[i, i] = intersection_number
            M[i, (i + 1) % n] = 1
            M[(i + 1) % n, i] = 1

        return M

    def compute_picard_generators(self):
        rays = self.rays

        D_description = []

        D_zero = [0, 0]
        D_one = [0, 0]

        for i in range(2, len(rays)):
            D_one.append(-int(np.inner([1, 0], rays[i])))
            D_zero.append(int(np.inner([rays[1][1], -1], rays[i])))

        D_description.append(D_zero)
        D_description.append(D_one)

        for j in range(2, len(rays) - 2):
            desc = [0] * len(rays)
            desc[j] = 1
            D_description.append(desc)

        return D_description


