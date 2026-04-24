# gui.py
import numpy as np
import string
import matplotlib.pyplot as plt
from Kaehler_Cone_Inequalities import (
    compute_nonredundant_inequalities,
    get_inequalities
)
from plotting import Plotter
from evaluator import evaluate_point
from matplotlib.widgets import Slider


class App:
    """
    Interactive controller for J(alpha, beta | geometry)
    """

    def __init__(self, resolution=100):

        self.geometry = None

        self.alpha = None
        self.beta = None

        self.alpha_sliders = {}
        self.beta_sliders = {}

        self.i = 0
        self.j = 0

        self.resolution = resolution
        self.axis_grid = np.linspace(np.finfo(float).eps, 1.0, resolution)

        self.plotter = Plotter()

    # ----------------------------
    # GEOMETRY SETUP
    # ----------------------------

    def set_geometry(self, geometry):
        self.geometry = geometry

        dim = geometry.divisor_dimension

        self.alpha = np.ones(dim)
        self.beta = np.ones(dim)

        self.i = 0
        self.j = 0

        self.create_sliders()

        raw_ineq = compute_nonredundant_inequalities(self.geometry)
        ineq_strings = get_inequalities(raw_ineq)
        self.plotter.set_inequalities(ineq_strings)
        self.update_plot()

    # ----------------------------
    # SLIDERS
    # ----------------------------

    def create_sliders(self):

        fig = self.plotter.fig

        # clear old sliders
        self.alpha_sliders.clear()
        self.beta_sliders.clear()

        y = 0.02
        h = 0.02
        spacing = 0.03

        dim = len(self.alpha)
        alphabet = string.ascii_lowercase
        # -------------------------
        # alpha sliders
        # -------------------------
        for k in range(dim):

            if k == self.i:
                continue

            ax_slider = fig.add_axes([0.05, y, 0.35, h])

            slider = Slider(
                ax=ax_slider,
                label=f"α[{alphabet[k]}]",
                valmin=1e-6,
                valmax=1.0,
                valinit=self.alpha[k]
            )

            slider.on_changed(self._make_alpha_callback(k))
            self.alpha_sliders[k] = slider

            y += spacing

        for k in range(dim):

            if k == self.j:
                continue

            ax_slider = fig.add_axes([0.05, y, 0.35, h])

            slider = Slider(
                ax=ax_slider,
                label=f"β[{alphabet[k]}]",
                valmin=1e-6,
                valmax=1.0,
                valinit=self.alpha[k]
            )

            slider.on_changed(self._make_beta_callback(k))
            self.beta_sliders[k] = slider

            y += spacing





        plt.draw()

    def _make_alpha_callback(self, k):
        def update(val):
            self.alpha[k] = max(val, 1e-6)
            self.update_plot()
        return update

    def _make_beta_callback(self, k):
        def update(val):
            self.beta[k] = max(val, 1e-6)
            self.update_plot()
        return update

    # ----------------------------
    # AXES SELECTION
    # ----------------------------

    def set_axes(self, i, j):
        if i == j:
            raise ValueError("i and j must be different")

        self.i = i
        self.j = j

        self.create_sliders()
        self.update_plot()

    # ----------------------------
    # COMPUTATION
    # ----------------------------


    def evaluate_slice(self):

        n = len(self.axis_grid)
        X, Y = np.meshgrid(self.axis_grid, self.axis_grid, indexing="ij")

        J = np.zeros((n, n))
        mask = np.zeros((n, n), dtype=bool)

        base_alpha = self.alpha.copy()
        base_beta = self.beta.copy()

        for a in range(n):
            for b in range(n):

                alpha = base_alpha.copy()
                beta = base_beta.copy()   # unchanged

                alpha[self.i] = X[a, b]
                beta[self.j] = Y[a, b]

                val, valid = evaluate_point(
                    (alpha, beta),
                    self.geometry
                )

                J[a, b] = val
                mask[a, b] = valid
                #print(f"J[a,b] = {J[a, b]}" + f"a = {a} : b = {b}")
                #print(f"alpha = {alpha}")
                #print(f"beta = {beta}")
        #print(f"J = {J}")
        return X, Y, J, mask
    # ----------------------------
    # UPDATE LOOP
    # ----------------------------

    def update_plot(self):

        if self.geometry is None:
            return

        X, Y, J, mask = self.evaluate_slice()

        if self.plotter.im is None:
            self.plotter.plot_values(
                X, Y, J,
                mask=mask,
                alpha_idx=self.i,
                beta_idx=self.j
            )
        else:
            self.plotter.update_figure(J, mask=mask)

        self.plotter.update_equations(
            self.alpha,
            self.beta,
            self.i,
            self.j
        )

    # ----------------------------
    # RUN
    # ----------------------------

    def run(self):
        plt.show()