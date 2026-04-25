# gui.py
import numpy as np
import string
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from Kaehler_Cone_Inequalities import (
    compute_nonredundant_inequalities,
    get_inequalities,
)
from plotting import Plotter
from evaluator import evaluate_point, EQUATIONS


class AxisChoice:
    def __init__(self, vector: str, k: int):
        assert vector in ("alpha", "beta")
        self.vector = vector
        self.k = k

    def __eq__(self, other):
        return self.vector == other.vector and self.k == other.k

    def label(self, alphabet):
        prefix = "α" if self.vector == "alpha" else "β"
        return f"{prefix}[{alphabet[self.k]}]"


class App:
    """
    Interactive controller for equation(alpha, beta | geometry).
    Either axis can be any α[k] or β[k].
    Equation can be switched via radio buttons.
    """

    def __init__(self, resolution=100):
        self.geometry = None

        self.alpha = None
        self.beta  = None

        self.alpha_sliders = {}
        self.beta_sliders  = {}

        self._x_radio  = None
        self._y_radio  = None
        self._eq_radio = None

        self.x_axis = AxisChoice("alpha", 0)
        self.y_axis = AxisChoice("beta",  0)

        self.equation_name = next(iter(EQUATIONS))

        self._pending_rebuild = False

        self.resolution = resolution
        self.axis_grid  = np.linspace(np.finfo(float).eps, 1.0, resolution)

        self.plotter = Plotter()
        self.plotter.on_resolution_change = self.set_resolution

    # -------------------------------------------------------------------
    # GEOMETRY SETUP
    # -------------------------------------------------------------------

    def set_geometry(self, geometry):
        self.geometry = geometry

        self.plotter.plot_fan(geometry.rays)

        dim = geometry.divisor_dimension

        self.alpha = np.ones(dim)
        self.beta  = np.ones(dim)

        self.x_axis = AxisChoice("alpha", 0)
        self.y_axis = AxisChoice("beta", 0) if dim > 1 else AxisChoice("alpha", 1)

        self._rebuild_controls()

        raw_ineq     = compute_nonredundant_inequalities(self.geometry)
        ineq_strings = get_inequalities(raw_ineq)
        self.plotter.set_inequalities(ineq_strings)
        self.update_plot()

    # -------------------------------------------------------------------
    # WIDGET TEARDOWN
    # -------------------------------------------------------------------

    def _remove_widget(self, widget):
        """
        Safely disconnect and remove a Slider or RadioButtons widget.
        Releases any mouse grab first to avoid 'Another Axes already
        grabs mouse input'.
        """
        if widget is None:
            return
        try:
            canvas = self.plotter.fig.canvas
            if hasattr(canvas, "mouse_grabber") and canvas.mouse_grabber is widget.ax:
                canvas.release_mouse(widget.ax)
        except Exception:
            pass
        try:
            widget.disconnect_events()
        except Exception:
            pass
        try:
            widget.ax.remove()
        except Exception:
            pass

    def _remove_all_widgets(self):
        self._remove_widget(self._x_radio)
        self._remove_widget(self._y_radio)
        self._remove_widget(self._eq_radio)
        for slider in list(self.alpha_sliders.values()) + list(self.beta_sliders.values()):
            self._remove_widget(slider)
        self.alpha_sliders.clear()
        self.beta_sliders.clear()
        self._x_radio  = None
        self._y_radio  = None
        self._eq_radio = None

    # -------------------------------------------------------------------
    # HELPERS
    # -------------------------------------------------------------------

    def _all_labels(self):
        dim      = len(self.alpha)
        alphabet = string.ascii_lowercase
        return (
                [f"α[{alphabet[k]}]" for k in range(dim)] +
                [f"β[{alphabet[k]}]" for k in range(dim)]
        )

    def _label_to_axis(self, label):
        dim      = len(self.alpha)
        alphabet = string.ascii_lowercase
        for k in range(dim):
            if label == f"α[{alphabet[k]}]":
                return AxisChoice("alpha", k)
            if label == f"β[{alphabet[k]}]":
                return AxisChoice("beta", k)
        raise ValueError(f"Unknown axis label: {label!r}")

    def _active_index(self, choice):
        dim = len(self.alpha)
        return choice.k if choice.vector == "alpha" else dim + choice.k

    def _axes_conflict(self, x: AxisChoice, y: AxisChoice):
        return x == y

    # -------------------------------------------------------------------
    # CONTROL REBUILD
    # -------------------------------------------------------------------

    def _rebuild_controls(self):
        self._remove_all_widgets()

        fig      = self.plotter.fig
        dim      = len(self.alpha)
        alphabet = string.ascii_lowercase
        labels   = self._all_labels()
        eq_names = list(EQUATIONS.keys())

        # -- Layout constants --------------------------------------------
        #
        #   [eq radio]  [x radio]  [y radio]      <- top strip
        #   [sliders, bottom-anchored]             <- below plot
        #
        n_axis_labels = 2 * dim
        n_eq_labels   = len(eq_names)

        r_width    = 0.05
        gap        = 0.02
        axis_h     = max(0.028 * n_axis_labels, 0.08)
        eq_h       = max(0.040 * n_eq_labels,   0.08)
        top        = 1.0
        left_eq    = 0.0
        left_x     = left_eq + r_width + gap + 0.15
        left_y     = left_x  + r_width + gap + 0.3
        slider_left = 0.05
        slider_w   = 0.30
        slider_h   = 0.025
        spacing    = 0.035

        # -- Equation radio ----------------------------------------------
        ax_eq = fig.add_axes([left_eq, top - eq_h, r_width, eq_h])
        ax_eq.set_title("equation", fontsize=8, pad=2)
        active_eq      = eq_names.index(self.equation_name) if self.equation_name in eq_names else 0
        self._eq_radio = RadioButtons(ax_eq, eq_names, active=active_eq)
        self._eq_radio.on_clicked(self._on_equation_changed)

        # -- X-axis radio ------------------------------------------------
        ax_x = fig.add_axes([left_x, top - axis_h, r_width, axis_h])
        ax_x.set_title("x-axis", fontsize=8, pad=2)
        self._x_radio = RadioButtons(ax_x, labels, active=self._active_index(self.x_axis))
        self._x_radio.on_clicked(self._on_x_axis_changed)

        # -- Y-axis radio ------------------------------------------------
        ax_y = fig.add_axes([left_y, top - axis_h, r_width, axis_h])
        ax_y.set_title("y-axis", fontsize=8, pad=2)
        self._y_radio = RadioButtons(ax_y, labels, active=self._active_index(self.y_axis))
        self._y_radio.on_clicked(self._on_y_axis_changed)

        # -- Sliders -----------------------------------------------------
        y = 0.2

        for k in range(dim):
            if AxisChoice("alpha", k) in (self.x_axis, self.y_axis):
                continue
            ax_s   = fig.add_axes([slider_left, y, slider_w, slider_h])
            slider = Slider(ax=ax_s, label=f"α[{alphabet[k]}]",
                            valmin=1e-6, valmax=1.0, valinit=self.alpha[k])
            slider.on_changed(self._make_alpha_callback(k))
            self.alpha_sliders[k] = slider
            y -= spacing

        for k in range(dim):
            if AxisChoice("beta", k) in (self.x_axis, self.y_axis):
                continue
            ax_s   = fig.add_axes([slider_left, y, slider_w, slider_h])
            slider = Slider(ax=ax_s, label=f"β[{alphabet[k]}]",
                            valmin=1e-6, valmax=1.0, valinit=self.beta[k])
            slider.on_changed(self._make_beta_callback(k))
            self.beta_sliders[k] = slider
            y -= spacing

        self._pending_rebuild = False
        fig.canvas.draw_idle()

    def _schedule_rebuild(self):
        """Defer _rebuild_controls to after the current event returns."""
        if self._pending_rebuild:
            return
        self._pending_rebuild = True
        fig = self.plotter.fig

        def _deferred(_):
            fig.canvas.mpl_disconnect(self._rebuild_cid)
            self._rebuild_controls()
            self.update_plot()

        self._rebuild_cid = fig.canvas.mpl_connect("draw_event", _deferred)
        fig.canvas.draw_idle()

    # -------------------------------------------------------------------
    # RADIO CALLBACKS
    # -------------------------------------------------------------------

    def set_resolution(self, res):
        self.resolution = res
        self.axis_grid = np.linspace(
            np.finfo(float).eps,
            1.0,
            self.resolution
        )
        self.update_plot()

    def _on_equation_changed(self, label):
        self.equation_name = label
        # No layout change needed — just recompute
        self.update_plot()

    def _on_x_axis_changed(self, label):
        new_x = self._label_to_axis(label)
        if self._axes_conflict(new_x, self.y_axis):
            self.y_axis = self.x_axis   # swap
        self.x_axis = new_x
        self._schedule_rebuild()

    def _on_y_axis_changed(self, label):
        new_y = self._label_to_axis(label)
        if self._axes_conflict(self.x_axis, new_y):
            self.x_axis = self.y_axis   # swap
        self.y_axis = new_y
        self._schedule_rebuild()

    # -------------------------------------------------------------------
    # SLIDER CALLBACKS
    # -------------------------------------------------------------------

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

    # -------------------------------------------------------------------
    # PROGRAMMATIC AXIS SELECTION
    # -------------------------------------------------------------------

    def set_axes(self, x_axis: AxisChoice, y_axis: AxisChoice):
        if self._axes_conflict(x_axis, y_axis):
            raise ValueError("x_axis and y_axis must refer to different coefficients")
        self.x_axis = x_axis
        self.y_axis = y_axis
        self._rebuild_controls()
        self.update_plot()

    # -------------------------------------------------------------------
    # COMPUTATION
    # -------------------------------------------------------------------

    def evaluate_slice(self):
        n          = len(self.axis_grid)
        X, Y       = np.meshgrid(self.axis_grid, self.axis_grid, indexing="ij")
        J          = np.zeros((n, n))
        mask       = np.zeros((n, n), dtype=bool)
        base_alpha = self.alpha.copy()
        base_beta  = self.beta.copy()

        for a in range(n):
            for b in range(n):
                alpha = base_alpha.copy()
                beta  = base_beta.copy()

                if self.x_axis.vector == "alpha":
                    alpha[self.x_axis.k] = X[a, b]
                else:
                    beta[self.x_axis.k]  = X[a, b]

                if self.y_axis.vector == "alpha":
                    alpha[self.y_axis.k] = Y[a, b]
                else:
                    beta[self.y_axis.k]  = Y[a, b]

                val, valid = evaluate_point(self.equation_name, (alpha, beta), self.geometry)
                J[a, b]    = val
                mask[a, b] = valid
        return X, Y, J, mask

    # -------------------------------------------------------------------
    # UPDATE LOOP
    # -------------------------------------------------------------------

    def update_plot(self):
        if self.geometry is None:
            return

        self.plotter.set_status("Computing...")
        self.plotter.fig.canvas.flush_events()

        X, Y, J, mask = self.evaluate_slice()

        self.plotter.set_status("Done")

        if self.plotter.im is None:
            self.plotter.plot_values(
                X, Y, J,
                mask=mask,
                x_axis=self.x_axis,
                y_axis=self.y_axis,
            )
        else:
            self.plotter.update_figure(J, mask=mask)

        self.plotter.update_equations(
            self.alpha,
            self.beta,
            self.x_axis,
            self.y_axis,
            self.equation_name,
        )

    # -------------------------------------------------------------------
    # RUN
    # -------------------------------------------------------------------

    def run(self):
        plt.show()