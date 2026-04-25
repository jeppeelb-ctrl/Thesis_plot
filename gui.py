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
from evaluator import evaluate_point


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
    Interactive controller for J(alpha, beta | geometry).
    Either axis can be any α[k] or β[k].
    """

    def __init__(self, resolution=100):

        self.geometry = None

        self.alpha = None
        self.beta  = None

        self.alpha_sliders = {}
        self.beta_sliders  = {}

        self._x_radio = None
        self._y_radio = None

        self.x_axis = AxisChoice("alpha", 0)
        self.y_axis = AxisChoice("beta",  0)

        self._pending_rebuild = False

        self.resolution = resolution
        self.axis_grid  = np.linspace(np.finfo(float).eps, 1.0, resolution)

        self.plotter = Plotter()

    # -------------------------------------------------------------------
    # GEOMETRY SETUP
    # -------------------------------------------------------------------

    def set_geometry(self, geometry):
        self.geometry = geometry

        dim = geometry.divisor_dimension

        self.alpha = np.ones(dim)
        self.beta  = np.ones(dim)

        self.x_axis = AxisChoice("alpha", 0)
        self.y_axis = AxisChoice("beta",  0) if dim > 0 else AxisChoice("alpha", min(1, dim - 1))

        self._rebuild_controls()

        raw_ineq    = compute_nonredundant_inequalities(self.geometry)
        ineq_strings = get_inequalities(raw_ineq)
        self.plotter.set_inequalities(ineq_strings)
        self.update_plot()

    # -------------------------------------------------------------------
    # WIDGET TEARDOWN HELPER
    # -------------------------------------------------------------------

    def _remove_widget(self, widget):
        """
        Safely disconnect and remove a Slider or RadioButtons widget.
        Releases any mouse grab before touching the axes so matplotlib
        does not throw 'Another Axes already grabs mouse input'.
        """
        if widget is None:
            return
        try:
            # Release mouse grab if this widget holds it
            canvas = self.plotter.fig.canvas
            if hasattr(canvas, "mouse_grabber") and canvas.mouse_grabber is widget.ax:
                canvas.release_mouse(widget.ax)
        except Exception:
            pass
        try:
            # Disconnect all callbacks registered by the widget
            widget.disconnect_events()
        except Exception:
            pass
        try:
            widget.ax.remove()
        except Exception:
            pass

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
        raise ValueError(f"Unknown axis label: {label}")

    def _active_index(self, choice):
        dim = len(self.alpha)
        return choice.k if choice.vector == "alpha" else dim + choice.k

    def _axes_conflict(self, x: AxisChoice, y: AxisChoice):
        return x == y

    # -------------------------------------------------------------------
    # CONTROL REBUILD
    # -------------------------------------------------------------------

    def _rebuild_controls(self):
        fig = self.plotter.fig

        # Safely tear down every existing widget before touching any axes
        self._remove_widget(self._x_radio)
        self._remove_widget(self._y_radio)

        for slider in list(self.alpha_sliders.values()) + list(self.beta_sliders.values()):
            self._remove_widget(slider)

        self.alpha_sliders.clear()
        self.beta_sliders.clear()
        self._x_radio = None
        self._y_radio = None

        dim      = len(self.alpha)
        n_labels = 2 * dim
        labels   = self._all_labels()
        alphabet = string.ascii_lowercase

        # Layout
        left     = 0
        r_width  = 0.13
        r_height = max(0.028 * n_labels, 0.08)
        gap      = 0.02
        slider_w = 0.30
        slider_h = 0.025
        spacing  = 0.035
        top      = 1

        # X-axis radio
        ax_x = fig.add_axes([left, top - r_height, r_width, r_height])
        ax_x.set_title("x-axis", fontsize=8, pad=2)
        self._x_radio = RadioButtons(ax_x, labels, active=self._active_index(self.x_axis))
        self._x_radio.on_clicked(self._on_x_axis_changed)

        # Y-axis radio
        ax_y = fig.add_axes([left + r_width + gap, top - r_height, r_width, r_height])
        ax_y.set_title("y-axis", fontsize=8, pad=2)
        self._y_radio = RadioButtons(ax_y, labels, active=self._active_index(self.y_axis))
        self._y_radio.on_clicked(self._on_y_axis_changed)

        # Sliders for every free coefficient
        #y = top - r_height - spacing
        slider_left = 0.05
        y = 0.15
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

    def _on_x_axis_changed(self, label):
        new_x = self._label_to_axis(label)
        if self._axes_conflict(new_x, self.y_axis):
            self.y_axis = self.x_axis
        self.x_axis = new_x
        self._schedule_rebuild()

    def _on_y_axis_changed(self, label):
        new_y = self._label_to_axis(label)
        if self._axes_conflict(self.x_axis, new_y):
            self.x_axis = self.y_axis
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
        n        = len(self.axis_grid)
        X, Y     = np.meshgrid(self.axis_grid, self.axis_grid, indexing="ij")
        J        = np.zeros((n, n))
        mask     = np.zeros((n, n), dtype=bool)
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

                val, valid    = evaluate_point((alpha, beta), self.geometry)
                J[a, b]       = val
                mask[a, b]    = valid

        return X, Y, J, mask

    # -------------------------------------------------------------------
    # UPDATE LOOP
    # -------------------------------------------------------------------

    def update_plot(self):
        if self.geometry is None:
            return

        X, Y, J, mask = self.evaluate_slice()

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
        )

    # -------------------------------------------------------------------
    # RUN
    # -------------------------------------------------------------------

    def run(self):
        plt.show()