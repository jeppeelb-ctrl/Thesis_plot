# geometry_manager.py
import numpy as np
from geometry import GeometryBuilder
from io_utils import parse_rays, get_surface_rays
from evaluator import (
    intersection_calculator,
    compute_anti_canonical,
)


class GeometryManager:

    SURFACE_GROUPS = {
        "P2": [
            ("P2Sigma0", "Cone 0:  ([0,1], [1,0])"),
            ("P2Sigma1", "Cone 1:  ([1,0], [-1,-1])"),
            ("P2Sigma2", "Cone 2:  ([-1,-1], [0,-1])"),
        ],
        "P1xP1": [
            ("P1xP1", "P1 × P1"),
        ],
        "Hirzebruch": [
            ("Hirzebruch_1", "r = 1"),
            ("Hirzebruch_2", "r = 2"),
            ("Hirzebruch_3", "r = 3"),
            ("Hirzebruch_4", "r = 4"),
        ],
    }

    # Rays that belong to the minimal base fan of each surface family.
    # These can never be blown down regardless of intersection numbers.
    _BASE_RAYS = {
        "P2":         {(0, 1), (1, 0), (-1, -1)},
        "P1xP1":      {(0, 1), (1, 0), (0, -1), (-1, 0)},
        "Hirzebruch": {(0, 1), (1, 0), (0, -1)},
        # (-1, r) is added dynamically in select_surface since r varies
    }

    def __init__(self):
        self.current     = None
        self._rays       = []
        self._base_rays  = set()   # populated in select_surface

    # ------------------------------------------------------------------
    # Surface selection
    # ------------------------------------------------------------------

    def select_surface(self, group_name, variant_key):
        if group_name == "Hirzebruch":
            r            = int(variant_key.split("_")[1])
            ray_str      = get_surface_rays("Hirzebruch")
            ray_str      = ray_str[:30] + str(r) + ray_str[31:]
            # base rays include the r-dependent ray
            self._base_rays = self._BASE_RAYS["Hirzebruch"] | {(-1, r)}
        else:
            self._base_rays = self._BASE_RAYS[group_name]

        self._rays   = parse_rays(ray_str if group_name == "Hirzebruch"
                                  else get_surface_rays(variant_key))
        self.current = GeometryBuilder(self._rays).build()

    # ------------------------------------------------------------------
    # Blow-up
    # ------------------------------------------------------------------

    def blow_up(self, cone_index):
        if self.current is None:
            raise RuntimeError("No surface selected.")

        rays = list(self._rays)
        n    = len(rays)
        i    = cone_index % n
        j    = (cone_index + 1) % n

        new_ray   = (
            rays[i][0] + rays[j][0],
            rays[i][1] + rays[j][1],
        )
        insert_at = j if j > i else i + 1
        rays.insert(insert_at, new_ray)

        self._rays   = rays
        self.current = GeometryBuilder(self._rays).build()

    # ------------------------------------------------------------------
    # Blow-down validity
    # ------------------------------------------------------------------

    def _is_base_ray(self, ray):
        """True if this ray is part of the minimal base fan."""
        return tuple(ray) in self._base_rays

    def _is_minus_one_curve(self, ray_index):
        """
        D_i is a (-1)-curve iff:
          1.  D_i^2 = -1
          2.  D_i · (-K_X) = 1
        """
        surface = self.current
        n       = len(self._rays)

        div              = [0] * n
        div[ray_index % n] = 1

        if intersection_calculator(div, div, surface) != -1:
            return False

        antican = compute_anti_canonical(surface)
        if intersection_calculator(div, antican, surface) != 1:
            return False

        return True

    def _can_blow_down(self, ray_index):
        """
        A ray can be blown down iff it is a (-1)-curve
        and does not belong to the base fan of the chosen surface family.
        """
        ray = self._rays[ray_index % len(self._rays)]
        return (
                not self._is_base_ray(ray)
                and self._is_minus_one_curve(ray_index)
        )

    def blowdown_valid_indices(self):
        return [i for i in range(len(self._rays)) if self._can_blow_down(i)]

    # ------------------------------------------------------------------
    # Blow-down
    # ------------------------------------------------------------------

    def blow_down(self, ray_index):
        if self.current is None:
            raise RuntimeError("No surface selected.")

        ray = self._rays[ray_index % len(self._rays)]

        if self._is_base_ray(ray):
            raise ValueError(
                f"Ray {ray} belongs to the base fan and cannot be blown down."
            )
        if not self._is_minus_one_curve(ray_index):
            raise ValueError(
                f"Ray {ray} is not a (−1)-curve and cannot be blown down."
            )

        rays = list(self._rays)
        rays.pop(ray_index % len(rays))
        self._rays   = rays
        self.current = GeometryBuilder(self._rays).build()

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def n_rays(self):
        return len(self._rays)

    @property
    def ray_labels(self):
        """
        Label each ray. Marks:
          (−1)       — valid blow-down target
          [base]     — belongs to the minimal fan, cannot be removed
        """
        valid = set(self.blowdown_valid_indices())
        labels = []
        for i, r in enumerate(self._rays):
            if i in valid:
                suffix = "  (−1)"
            elif self._is_base_ray(r):
                suffix = "  [base]"
            else:
                suffix = ""
            labels.append(f"Ray {i}:  {r}{suffix}")
        return labels

    @property
    def cone_labels(self):
        n = len(self._rays)
        return [
            f"Cone {i}:  {self._rays[i]} ∧ {self._rays[(i + 1) % n]}"
            for i in range(n)
        ]

    @property
    def intersection_matrix(self):
        return self.current.intersection_matrix if self.current is not None else None