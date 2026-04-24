import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import TextBox
import string


class Plotter:

    def __init__(self):

        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)

        self.im = None
        self.cbar = None

        # colormap
        self.cmap = plt.cm.RdYlGn.copy()
        self.cmap.set_bad(color="grey")

        # normalization
        self.vmin = -10
        self.vmax = 10
        self.norm = colors.TwoSlopeNorm(vmin=self.vmin, vcenter=0, vmax=self.vmax)

        # -------------------------
        # TEXTBOXES (vmin/vmax)
        # -------------------------
        ax_vmin = self.fig.add_axes([0.9, 0.5, 0.05, 0.05])
        ax_vmax = self.fig.add_axes([0.9, 0.3, 0.05, 0.05])

        self.tb_vmin = TextBox(ax_vmin, "vmin", initial=str(self.vmin))
        self.tb_vmax = TextBox(ax_vmax, "vmax", initial=str(self.vmax))

        self.tb_vmin.on_submit(self._update_norm_from_text)
        self.tb_vmax.on_submit(self._update_norm_from_text)

        # -------------------------
        # α and β coefficients
        # -------------------------
        self.alpha_coeffs = [1, 2, 3]
        self.beta_coeffs = [3, 1, 4]

        self.alpha_idx = 0
        self.beta_idx = 0

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
    # INITIAL PLOT
    # =========================================================
    def plot_values(self, X, Y, J, mask=None, alpha_idx=0, beta_idx=0):

        self.alpha_idx = alpha_idx
        self.beta_idx = beta_idx

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

        alphabet = string.ascii_lowercase

        self.ax.set_title(rf"J(α,β)")
        self.ax.set_xlabel(rf"β{alphabet[alpha_idx]}")
        self.ax.set_ylabel(rf"α{alphabet[beta_idx]}")

        self.update_equations(
            self.alpha_coeffs,
            self.beta_coeffs,
            self.alpha_idx,
            self.beta_idx
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
            self.alpha_idx,
            self.beta_idx
        )

        self.fig.canvas.draw_idle()

    # =========================================================
    # EQUATIONS
    # =========================================================
    def update_equations(self, alpha, beta, alpha_idx, beta_idx):

        def build(coeffs, active_idx, symbol):

            terms = []

            for i, a in enumerate(coeffs):

                if i == active_idx:
                    terms.append(rf"{symbol}_{{{i}}}")
                else:
                    terms.append(rf"{a:.3f} {symbol}_{{{i}}}")

            return " + ".join(terms)

        alpha_expr = build(alpha, alpha_idx, r"g")
        beta_expr = build(beta, beta_idx, r"g")

        self.alpha_text.set_text(rf"$\alpha = {alpha_expr}$")
        self.beta_text.set_text(rf"$\beta = {beta_expr}$")

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

    # =========================================================
    # ARGMIN
    # =========================================================
    def plot_argmin(self, X, Y, J, mask=None):

        if mask is not None:
            J = np.ma.array(J, mask=~mask)
            idx = np.unravel_index(np.argmin(J), J.shape)
        else:
            idx = np.unravel_index(np.argmin(J), J.shape)

        x_min = X[idx]
        y_min = Y[idx]

        self.ax.plot(x_min, y_min, "ro")
        self.fig.canvas.draw_idle()