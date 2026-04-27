# plotting.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import TextBox, Button
import string


class Plotter:

    def __init__(self):

        self.equation_name = None
        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)

        self.im   = None
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
        # -------------------------
        self.ax_fan = self.fig.add_axes([0.6, 0.1, 0.1, 0.1])
        self.ax_fan.set_aspect("equal")
        self.ax_fan.axis("off")

        # -------------------------
        # TEXTBOXES (vmin/vmax, resolution)
        # -------------------------
        ax_vmin = self.fig.add_axes([0.9, 0.5, 0.05, 0.05])
        ax_vmax = self.fig.add_axes([0.9, 0.3, 0.05, 0.05])
        ax_res  = self.fig.add_axes([0.9, 0.1, 0.05, 0.05])

        self.tb_vmin = TextBox(ax_vmin, "vmin",       initial=str(self.vmin))
        self.tb_vmax = TextBox(ax_vmax, "vmax",       initial=str(self.vmax))
        self.tb_res  = TextBox(ax_res,  "resolution", initial=str(self.resolution))

        self.tb_vmin.on_submit(self._update_norm_from_text)
        self.tb_vmax.on_submit(self._update_norm_from_text)
        self.tb_res.on_submit(self._update_resolution)

        # -------------------------
        # INEQUALITY TEXT
        # -------------------------

        self.ineq_text = None
        self._ineq_texts = []

        # -------------------------
        # DIVISOR OVERLAY
        # -------------------------
        self._show_divisors  = True
        self._divisor_im     = None   # imshow overlay
        self._divisor_labels = []     # text artists
        self._div_X_cache    = None
        self._div_Y_cache    = None
        self._div_grid_cache = None
        self._div_mask_cache = None

        ax_toggle = self.fig.add_axes([0.45, 0.12, 0.1, 0.04])
        self._div_toggle_btn = Button(ax_toggle, "Divisors: ON")
        self._div_toggle_btn.on_clicked(self._on_divisor_toggle)

        # -------------------------
        # ZERO LOCUS OVERLAY
        # -------------------------
        self._show_zero_locus    = True
        self._zero_locus_contour = None   # single ContourSet
        self._zl_X_cache         = None
        self._zl_Y_cache         = None
        self._zl_J_cache         = None
        self._zl_mask_cache      = None

        ax_zl_toggle = self.fig.add_axes([0.45, 0.17, 0.1, 0.04])
        self._zl_toggle_btn = Button(ax_zl_toggle, "Zero locus: ON")
        self._zl_toggle_btn.on_clicked(self._on_zero_locus_toggle)

        # -------------------------
        # STATUS TEXT
        # -------------------------
        self.status_text = self.ax.text(
            0.98, 0.02, "",
            transform=self.ax.transAxes,
            ha="right", va="bottom",
            fontsize=12,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
        )

        # -------------------------
        # α and β coefficients
        # -------------------------
        self.alpha_coeffs = [1, 2, 3]
        self.beta_coeffs  = [3, 1, 4]

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
        alphabet = string.ascii_lowercase
        prefix   = "α" if axis.vector == "alpha" else "β"
        return rf"${prefix}_{{{alphabet[axis.k]}}}$"

    # =========================================================
    # INITIAL PLOT
    # =========================================================

    def plot_values(self, X, Y, J, mask=None, x_axis=None, y_axis=None):

        self.x_axis = x_axis
        self.y_axis = y_axis

        self._clear_divisor_overlay()
        self._clear_zero_locus()

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

        self._clear_divisor_overlay()
        self._clear_zero_locus()

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

    def _update_axis_labels(self):
        if self.x_axis is None or self.y_axis is None:
            return
        self.ax.set_xlabel(self._axis_label(self.y_axis))
        self.ax.set_ylabel(self._axis_label(self.x_axis))

    # =========================================================
    # EQUATIONS
    # =========================================================

    def update_equations(self, alpha, beta, x_axis, y_axis, equation_name=""):

        self.alpha_coeffs  = list(alpha)
        self.beta_coeffs   = list(beta)
        self.x_axis        = x_axis
        self.y_axis        = y_axis
        self.equation_name = equation_name

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
                    terms.append(rf"{symbol}_{{{i}}}")
                else:
                    terms.append(rf"{a:.3f}\,{symbol}_{{{i}}}")
            return " + ".join(terms)

        alpha_expr = build(alpha, active_alpha, r"g")
        beta_expr  = build(beta,  active_beta,  r"g")

        self.alpha_text.set_text(rf"$\alpha = {alpha_expr}$")
        self.beta_text.set_text( rf"$\beta  = {beta_expr}$")

        self.ax.set_title(equation_name)
        self._update_axis_labels()

    # =========================================================
    # INEQUALITIES (INJECTED FROM APP)
    # =========================================================

    def set_inequalities(self, ineq_strings):
        if not ineq_strings:
            if self.ineq_text is not None:
                self.ineq_text.remove()
                self.ineq_text = None
            for t in self._ineq_texts:
                try:
                    t.remove()
                except Exception:
                    pass
            x_pos = 1.3
            y_start = 1.1
            y_step  = 0.08

            t = self.ax.text(
                x_pos, y_start,
                rf"$\;No\;inequalities\;recieved$",
                transform=self.ax.transAxes,
                va="top", ha="right",
                fontsize=11,
                color="black",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1),
            )

            self._ineq_texts.append(t)
            self.fig.canvas.draw_idle()
            return

        # Remove old
        if self.ineq_text is not None:
            self.ineq_text.remove()
            self.ineq_text = None
        for t in self._ineq_texts:
            try:
                t.remove()
            except Exception:
                pass
        self._ineq_texts = []

        self._ineq_strings = ineq_strings   # store for recoloring

        # Place one Text per inequality so we can color them individually
        x_pos = 1.3
        y_start = 1.1
        y_step  = 0.08

        for i, s in enumerate(ineq_strings):
            t = self.ax.text(
                x_pos, y_start - i * y_step,
                rf"$\;{s}$",
                transform=self.ax.transAxes,
                va="top", ha="right",
                fontsize=11,
                color="black",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1),
                       )
            self._ineq_texts.append(t)

        self.fig.canvas.draw_idle()

    def update_inequality_colors(self, satisfied: list[bool]):
        """
        Color each inequality label green if satisfied, red if violated.
        satisfied[i] corresponds to self._ineq_texts[i].
        """
        for t, ok in zip(self._ineq_texts, satisfied):
            t.set_color("#2a7a2a" if ok else "#cc2222")
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
    # TEXTBOX CALLBACKS
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
        self.norm = colors.TwoSlopeNorm(vmin=self.vmin, vcenter=0, vmax=self.vmax)

        if self.im is not None:
            self.im.set_norm(self.norm)

        self.fig.canvas.draw_idle()

    def _update_resolution(self, text):
        try:
            res = max(5, int(text))
        except ValueError:
            return
        if hasattr(self, "on_resolution_change"):
            self.on_resolution_change(res)

    def set_status(self, status: str):
        self.status_text.set_text(status)
        self.fig.canvas.draw_idle()

    # =========================================================
    # DIVISOR TOGGLE
    # =========================================================

    def _on_divisor_toggle(self, _):
        self._show_divisors = not self._show_divisors

        if self._show_divisors:
            self._div_toggle_btn.label.set_text("Divisors: ON")
            if self._div_grid_cache is not None:
                self._draw_divisor_regions(
                    self._div_X_cache,
                    self._div_Y_cache,
                    self._div_grid_cache,
                    self._div_mask_cache,
                )
        else:
            self._div_toggle_btn.label.set_text("Divisors: OFF")
            self._clear_divisor_overlay()

        self.fig.canvas.draw_idle()

    # =========================================================
    # DIVISOR OVERLAY — public entry point
    # =========================================================

    def plot_divisor_regions(self, X, Y, div_grid, mask=None):
        """Cache the latest data and draw if toggled on."""
        self._div_X_cache    = X
        self._div_Y_cache    = Y
        self._div_grid_cache = div_grid
        self._div_mask_cache = mask

        if self._show_divisors:
            self._draw_divisor_regions(X, Y, div_grid, mask)

    # =========================================================
    # DIVISOR OVERLAY — internal clear / draw
    # =========================================================

    def _clear_divisor_overlay(self):
        if self._divisor_im is not None:
            try:
                self._divisor_im.remove()
            except Exception:
                pass
            self._divisor_im = None

        for txt in self._divisor_labels:
            try:
                txt.remove()
            except Exception:
                pass
        self._divisor_labels = []

    def _draw_divisor_regions(self, X, Y, div_grid, mask=None):
        self._clear_divisor_overlay()

        masked_grid = np.where(mask, div_grid, -1) if mask is not None else div_grid.copy()

        valid_divs = masked_grid[masked_grid >= 0]
        if valid_divs.size == 0:
            return

        unique_divs = sorted(set(valid_divs.ravel()))
        if len(unique_divs) <= 1:
            return

        # Semi-transparent imshow overlay — same extent/origin as the main plot.
        # div_grid is indexed [x, y] so .T gives imshow's expected [row=y, col=x].
        float_grid          = masked_grid.astype(float)
        float_grid[masked_grid < 0] = np.nan
        display_grid        = np.ma.array(float_grid, mask=(masked_grid < 0))

        cmap_div = plt.cm.get_cmap("tab10", len(unique_divs))

        self._divisor_im = self.ax.imshow(
            display_grid,
            origin="lower",
            extent=[X.min(), X.max(), Y.min(), Y.max()],
            aspect="auto",
            interpolation="nearest",
            cmap=cmap_div,
            alpha=0.25,
            vmin=unique_divs[0]  - 0.5,
            vmax=unique_divs[-1] + 0.5,
            zorder=2,
        )

        # Label each region at its centroid
        for d in unique_divs:
            # Only consider points that are both this divisor AND valid
            if mask is not None:
                region = (masked_grid == d) & mask
            else:
                region = (masked_grid == d)

            if not region.any():
                continue

            if region.sum() < 20:
                continue

            cy = X[region].mean()
            cx = Y[region].mean()

            txt = self.ax.text(
                cx, cy,
                rf"$D_{{{d}}}$",
                ha="center", va="center",
                fontsize=10, color="white",
                fontweight="bold",
                zorder=4,
            )
            self._divisor_labels.append(txt)


    # =========================================================
    # ZERO LOCUS TOGGLE
    # =========================================================

    def _on_zero_locus_toggle(self, _):
        self._show_zero_locus = not self._show_zero_locus

        if self._show_zero_locus:
            self._zl_toggle_btn.label.set_text("Zero locus: ON")
            if self._zl_J_cache is not None:
                self._draw_zero_locus(
                    self._zl_X_cache,
                    self._zl_Y_cache,
                    self._zl_J_cache,
                    self._zl_mask_cache,
                )
        else:
            self._zl_toggle_btn.label.set_text("Zero locus: OFF")
            self._clear_zero_locus()

        self.fig.canvas.draw_idle()

    # =========================================================
    # ZERO LOCUS OVERLAY
    # =========================================================

    def plot_zero_locus(self, X, Y, J, mask=None):
        """Cache and draw the zero locus {J = 0} if toggled on."""
        self._zl_X_cache    = X
        self._zl_Y_cache    = Y
        self._zl_J_cache    = J
        self._zl_mask_cache = mask

        if self._show_zero_locus:
            self._draw_zero_locus(X, Y, J, mask)

    def _clear_zero_locus(self):
        if self._zero_locus_contour is not None:
            try:
                self._zero_locus_contour.remove()
            except Exception:
                pass
            self._zero_locus_contour = None

    def _draw_zero_locus(self, X, Y, J, mask=None):
        self._clear_zero_locus()

        # Work on a copy so we don't mutate the cache
        J_plot = J.copy().astype(float)

        if mask is not None:
            J_plot[~mask] = np.nan

        # Check there are both positive and negative valid values —
        # if not, the zero locus doesn't cross this slice at all
        valid = J_plot[~np.isnan(J_plot)]
        if valid.size == 0 or valid.min() >= 0 or valid.max() <= 0:
            return

        y_vals = X[:, 0]
        x_vals = Y[0, :]

        try:
            # Single contour at level 0, drawn as a solid black line
            self._zero_locus_contour = self.ax.contour(
                x_vals, y_vals,
                J_plot,           # contour expects (ny, nx)
                levels=[0],
                colors="black",
                linewidths=2.0,
                linestyles="-",
                zorder=5,
            )
        except Exception as e:
            print(f"Zero locus contour failed: {e}")

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