# plotting.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import TextBox
import string


class Plotter:

    def __init__(self):

        self.equation_name = None
        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)

        self.im = None
        self.cbar = None

        # colormap
        self.cmap = plt.cm.RdYlGn.copy()
        self.cmap.set_bad(color="grey")

        # normalization
        self.vmin = -5
        self.vmax = 5
        self.norm = colors.TwoSlopeNorm(vmin=self.vmin, vcenter=0, vmax=self.vmax)

        self.resolution = 100

        # -------------------------
        # FAN PICTURE
        # ------------------------

        # side panel for fan
        self.ax_fan = self.fig.add_axes([0.6, 0.1 , 0.1, 0.1])  # [left, bottom, width, height]
        self.ax_fan.set_aspect("equal")
        self.ax_fan.axis("off")

        # -------------------------
        # TEXTBOXES (vmin/vmax), resolution and computing or not
        # -------------------------
        ax_vmin = self.fig.add_axes([0.9, 0.5, 0.05, 0.05])
        ax_vmax = self.fig.add_axes([0.9, 0.3, 0.05, 0.05])
        ax_res = self.fig.add_axes([0.9, 0.1, 0.05, 0.05])

        self.tb_vmin = TextBox(ax_vmin, "vmin", initial=str(self.vmin))
        self.tb_vmax = TextBox(ax_vmax, "vmax", initial=str(self.vmax))
        self.tb_res = TextBox(ax_res, "resolution", initial=str(self.resolution))

        self.tb_vmin.on_submit(self._update_norm_from_text)
        self.tb_vmax.on_submit(self._update_norm_from_text)
        self.tb_res.on_submit(self._update_resolution)

        self.status_text = self.ax.text(
            0.98, 0.02,
            "",
            transform=self.ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=12,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
        )

        # -------------------------
        # α and β coefficients
        # -------------------------
        self.alpha_coeffs = [1, 2, 3]
        self.beta_coeffs  = [3, 1, 4]

        # AxisChoice objects — set properly on first plot_values call
        self.x_axis = None
        self.y_axis = None

        self.ineq_text = None

        # -------------------------
        # TEXT DISPLAY
        # -------------------------
        self.alpha_text = self.ax.text(
            0.02, 0.98, "",
            transform=self.ax.transAxes,
            va="top", ha="left",
            fontsize=11
        )

        self.beta_text = self.ax.text(
            0.02, 0.88, "",
            transform=self.ax.transAxes,
            va="top", ha="left",
            fontsize=11
        )

    # =========================================================
    # HELPERS
    # =========================================================

    @staticmethod
    def _axis_label(axis):
        """Human-readable label for an AxisChoice, e.g. 'α[a]' or 'β[b]'."""
        alphabet = string.ascii_lowercase
        prefix = "α" if axis.vector == "alpha" else "β"
        return rf"${prefix}_{{{alphabet[axis.k]}}}$"

    # =========================================================
    # INITIAL PLOT
    # =========================================================

    def plot_values(self, X, Y, J, mask=None, x_axis=None, y_axis=None):

        self.x_axis = x_axis
        self.y_axis = y_axis

        if mask is not None:
            J = np.ma.array(J, mask=~mask)

        self.im = self.ax.imshow(
            J,
            origin="lower",
            extent=[X.min(), X.max(), Y.min(), Y.max()],
            aspect="auto",
            interpolation="nearest",
            cmap=self.cmap,
            norm=self.norm
        )

        if self.cbar is None:
            self.cbar = self.fig.colorbar(self.im, ax=self.ax)

        self._update_axis_labels()

        self.update_equations(
            self.alpha_coeffs,
            self.beta_coeffs,
            self.x_axis,
            self.y_axis,
            self.equation_name
        )

        plt.draw()

    # =========================================================
    # LIVE UPDATE
    # =========================================================

    def update_figure(self, J, mask=None):

        if self.im is None:
            raise RuntimeError("Plot not initialized.")

        if mask is not None:
            J = np.ma.array(J, mask=~mask)

        self.im.set_array(J)

        self.update_equations(
            self.alpha_coeffs,
            self.beta_coeffs,
            self.x_axis,
            self.y_axis,
            self.equation_name
        )

        self.fig.canvas.draw_idle()

    # =========================================================
    # AXIS LABELS
    # =========================================================

    def _update_axis_labels(self, equation_name=""):
        """Sync the plot title and x/y axis labels with the current AxisChoice pair."""
        if self.x_axis is None or self.y_axis is None:
            return

        x_lbl = self._axis_label(self.x_axis)
        y_lbl = self._axis_label(self.y_axis)

        self.ax.set_xlabel(x_lbl)
        self.ax.set_ylabel(y_lbl)

    # =========================================================
    # EQUATIONS
    # =========================================================

    def update_equations(self, alpha, beta, x_axis, y_axis, equation_name=""):
        """
        Display the current α and β vectors, highlighting whichever
        components are being used as plot axes.
        """
        self.alpha_coeffs = list(alpha)
        self.beta_coeffs  = list(beta)
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.equation_name = equation_name

        # Collect which α and β indices are plot axes (may be 0, 1, or 2 each)
        active_alpha = set()
        active_beta  = set()
        for ax_choice in (x_axis, y_axis):
            if ax_choice is None:
                continue
            if ax_choice.vector == "alpha":
                active_alpha.add(ax_choice.k)
            else:
                active_beta.add(ax_choice.k)

        def build(coeffs, active_set, symbol):
            terms = []
            for i, a in enumerate(coeffs):
                if i in active_set:
                    # highlighted: show bare variable, no fixed value
                    terms.append(rf"{symbol}_{{{i}}}")
                else:
                    terms.append(rf"{a:.3f}\,{symbol}_{{{i}}}")
            return " + ".join(terms)

        alphabet = string.ascii_lowercase
        alpha_expr = build(alpha, active_alpha, r"g")
        beta_expr  = build(beta,  active_beta,  r"g")

        self.alpha_text.set_text(rf"$\alpha = {alpha_expr}$")
        self.beta_text.set_text( rf"$\beta  = {beta_expr}$")

        # Keep axis labels in sync whenever equations are refreshed
        self.ax.set_title(equation_name)
        self._update_axis_labels()

    # =========================================================
    # INEQUALITIES (INJECTED FROM APP)
    # =========================================================

    def set_inequalities(self, ineq_strings):
        if not ineq_strings:
            print("No inequalities received!")
            return

        if self.ineq_text is not None:
            self.ineq_text.remove()

        text = "\n".join([f"$\\;{s}$" for s in ineq_strings])

        self.ineq_text = self.ax.text(
            1.3, 1.1,
            text,
            transform=self.ax.transAxes,
            va="top",
            ha="right",
            fontsize=11,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
        )

        self.fig.canvas.draw_idle()


    # =========================================================
    # FAN PICTURE
    # =========================================================

    def plot_fan(self, rays):

        self.ax_fan.clear()
        self.ax_fan.set_aspect("equal")
        self.ax_fan.axis("off")

        for x, y in rays:
            self.ax_fan.plot([0, x], [0, y], "k-")
            self.ax_fan.plot(x, y, "ko", markersize=3)

        self.ax_fan.set_xlim(-1.2, 1.2)
        self.ax_fan.set_ylim(-1.2, 1.2)

        self.fig.canvas.draw_idle()

    # =========================================================
    # TEXTBOX CALLBACK
    # =========================================================

    def _update_norm_from_text(self, _):

        try:
            vmin = int(self.tb_vmin.text)
            vmax = int(self.tb_vmax.text)
        except ValueError:
            return

        if vmin >= vmax:
            return

        self.vmin = vmin
        self.vmax = vmax

        self.norm = colors.TwoSlopeNorm(
            vmin=self.vmin,
            vcenter=0,
            vmax=self.vmax
        )

        if self.im is not None:
            self.im.set_norm(self.norm)

        self.fig.canvas.draw_idle()


    def _update_resolution(self, text):
        try:
            res = int(text)
            res = max(5, res)
        except ValueError:
            return

        if hasattr(self, "on_resolution_change"):
            self.on_resolution_change(res)


    def set_status(self, status: str):
        self.status_text.set_text(status)
        self.fig.canvas.draw_idle()

    # =========================================================
    # ARGMIN
    # =========================================================

    def plot_argmin(self, X, Y, J, mask=None):

        if mask is not None:
            J = np.ma.array(J, mask=~mask)

        idx   = np.unravel_index(np.argmin(J), J.shape)
        x_min = X[idx]
        y_min = Y[idx]

        self.ax.plot(x_min, y_min, "ro")
        self.fig.canvas.draw_idle()