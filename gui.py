# gui.py
import numpy as np
import string
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, Button
from Kaehler_Cone_Inequalities import (
    compute_nonredundant_inequalities,
    get_inequalities,
)
from evaluator import evaluate_point, EQUATIONS, is_kaehler
from geometry_manager import GeometryManager


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

    def __init__(self, resolution=100):

        self.geo_manager = GeometryManager()
        self.plotter     = None

        self.alpha = None
        self.beta  = None

        # Window 1: equation + axis radios
        self._controls_fig = None
        self._x_radio      = None
        self._y_radio      = None
        self._eq_radio     = None

        # Window 2: sliders + inequalities
        self._sliders_fig  = None
        self.alpha_sliders = {}
        self.beta_sliders  = {}
        self._ineq_strings = []
        self._ineq_texts   = []

        # Window 3: geometry (blow-up/down) + intersection matrix
        self._ctrl_fig     = None
        self._cone_radio   = None
        self._ray_radio    = None
        self._blowup_btn   = None
        self._blowdown_btn = None
        self._ctrl_status  = None
        self._matrix_fig   = None

        self.x_axis = AxisChoice("alpha", 0)
        self.y_axis = AxisChoice("beta",  0)

        self.equation_name  = next(iter(EQUATIONS))
        self._selected_cone = 0
        self._selected_ray  = 0

        self._pending_rebuild = False

        self.resolution = resolution
        self.axis_grid  = np.linspace(np.finfo(float).eps, 1.0, resolution)

    # -------------------------------------------------------------------
    # STARTUP DIALOGS
    # -------------------------------------------------------------------

    def _pick_surface(self):
        groups         = list(GeometryManager.SURFACE_GROUPS.keys())
        chosen_group   = [groups[0]]
        chosen_variant = [None]

        fig1 = plt.figure(figsize=(3.5, 2.5))
        fig1.suptitle("Select surface family", fontsize=11)
        ax_r1  = fig1.add_axes([0.1, 0.30, 0.8, 0.55])
        ax_b1  = fig1.add_axes([0.25, 0.07, 0.5, 0.16])
        radio1 = RadioButtons(ax_r1, groups, active=0)
        btn1   = Button(ax_b1, "Next →")

        def _confirm_group(_):
            chosen_group[0] = radio1.value_selected
            plt.close(fig1)

        btn1.on_clicked(_confirm_group)
        plt.show(block=True)

        variants       = GeometryManager.SURFACE_GROUPS[chosen_group[0]]
        variant_labels = [v[1] for v in variants]
        variant_keys   = [v[0] for v in variants]

        fig2 = plt.figure(figsize=(4.5, max(2.5, 0.5 + 0.4 * len(variants))))
        fig2.suptitle(f"Select {chosen_group[0]} variant", fontsize=11)
        ax_r2  = fig2.add_axes([0.1, 0.25, 0.8, 0.60])
        ax_b2  = fig2.add_axes([0.25, 0.05, 0.5, 0.14])
        radio2 = RadioButtons(ax_r2, variant_labels, active=0)
        btn2   = Button(ax_b2, "Confirm")

        def _confirm_variant(_):
            idx               = variant_labels.index(radio2.value_selected)
            chosen_variant[0] = variant_keys[idx]
            plt.close(fig2)

        btn2.on_clicked(_confirm_variant)
        plt.show(block=True)

        return chosen_group[0], chosen_variant[0]

    # -------------------------------------------------------------------
    # WINDOW 1: CONTROLS  (equation + axis radios)
    # -------------------------------------------------------------------

    def _build_controls_figure(self):
        if self._controls_fig is not None:
            try:
                plt.close(self._controls_fig)
            except Exception:
                pass

        dim      = len(self.alpha)
        labels   = self._all_labels()
        eq_names = list(EQUATIONS.keys())
        n_axis   = 2 * dim

        eq_h   = max(0.06 * len(eq_names), 0.10)
        axis_h = max(0.032 * n_axis,       0.15)
        fig_h  = max(3.0, eq_h * 8 + axis_h * 8 + 1.0)

        self._controls_fig = plt.figure("Controls", figsize=(4.5, fig_h))
        self._controls_fig.subplots_adjust(
            left=0.05, right=0.95, top=0.97, bottom=0.03
        )

        r_w  = 0.42
        gap  = 0.06
        top  = 0.97

        # Equation radio
        ax_eq = self._controls_fig.add_axes([0.05, top - eq_h, 0.88, eq_h])
        ax_eq.set_title("Equation", fontsize=9, pad=2)
        active_eq      = eq_names.index(self.equation_name) if self.equation_name in eq_names else 0
        self._eq_radio = RadioButtons(ax_eq, eq_names, active=active_eq)
        self._eq_radio.on_clicked(self._on_equation_changed)

        # X-axis radio
        axis_top = top - eq_h - 0.04
        ax_x = self._controls_fig.add_axes([0.05, axis_top - axis_h, r_w, axis_h])
        ax_x.set_title("x-axis", fontsize=9, pad=2)
        self._x_radio = RadioButtons(ax_x, labels, active=self._active_index(self.x_axis))
        self._x_radio.on_clicked(self._on_x_axis_changed)

        # Y-axis radio
        ax_y = self._controls_fig.add_axes(
            [0.05 + r_w + gap, axis_top - axis_h, r_w, axis_h]
        )
        ax_y.set_title("y-axis", fontsize=9, pad=2)
        self._y_radio = RadioButtons(ax_y, labels, active=self._active_index(self.y_axis))
        self._y_radio.on_clicked(self._on_y_axis_changed)

        self._controls_fig.show()

    # -------------------------------------------------------------------
    # WINDOW 2: SLIDERS + INEQUALITIES
    # -------------------------------------------------------------------

    def _build_sliders_figure(self):
        if self._sliders_fig is not None:
            try:
                plt.close(self._sliders_fig)
            except Exception:
                pass

        dim      = len(self.alpha)
        alphabet = string.ascii_lowercase
        n_sliders = 2 * dim - 2
        n_ineq    = len(self._ineq_strings)

        fig_h = max(3.5,
                    0.40 * n_sliders
                    + 0.30 * n_ineq
                    + 1.2)

        self._sliders_fig = plt.figure("Sliders & conditions", figsize=(5.0, fig_h))
        self._sliders_fig.subplots_adjust(
            left=0.05, right=0.95, top=0.97, bottom=0.03
        )

        slider_h    = 0.045
        spacing     = max(0.055, 0.85 / max(n_sliders, 1))
        slider_left = 0.12
        slider_w    = 0.80
        top         = 0.95

        y = top
        self.alpha_sliders.clear()
        self.beta_sliders.clear()

        for k in range(dim):
            if AxisChoice("alpha", k) in (self.x_axis, self.y_axis):
                continue
            ax_s   = self._sliders_fig.add_axes([slider_left, y - slider_h, slider_w, slider_h])
            slider = Slider(ax=ax_s, label=f"α[{alphabet[k]}]",
                            valmin=1e-6, valmax=1.0, valinit=self.alpha[k])
            slider.on_changed(self._make_alpha_callback(k))
            self.alpha_sliders[k] = slider
            y -= spacing

        for k in range(dim):
            if AxisChoice("beta", k) in (self.x_axis, self.y_axis):
                continue
            ax_s   = self._sliders_fig.add_axes([slider_left, y - slider_h, slider_w, slider_h])
            slider = Slider(ax=ax_s, label=f"β[{alphabet[k]}]",
                            valmin=1e-6, valmax=1.0, valinit=self.beta[k])
            slider.on_changed(self._make_beta_callback(k))
            self.beta_sliders[k] = slider
            y -= spacing

        # Inequalities below sliders
        self._ineq_texts = []

        if self._ineq_strings:
            ineq_top  = y - 0.04
            ineq_step = min(0.07, (ineq_top - 0.02) / max(len(self._ineq_strings), 1))

            self._sliders_fig.text(
                0.5, ineq_top,
                "Kähler cone conditions",
                ha="center", va="top",
                fontsize=9, fontstyle="italic",
                transform=self._sliders_fig.transFigure,
            )

            for i, s in enumerate(self._ineq_strings):
                t = self._sliders_fig.text(
                    0.5,
                    ineq_top - 0.05 - i * ineq_step,
                    rf"${s}$",
                    ha="center", va="top",
                    fontsize=11,
                    color="black",
                    transform=self._sliders_fig.transFigure,
                    )
                self._ineq_texts.append(t)

        self._sliders_fig.show()

    def _remove_slider_widgets(self):
        for slider in list(self.alpha_sliders.values()) + list(self.beta_sliders.values()):
            if slider is None:
                continue
            try:
                if self._sliders_fig is not None:
                    canvas = self._sliders_fig.canvas
                    if hasattr(canvas, "mouse_grabber") and canvas.mouse_grabber is slider.ax:
                        canvas.release_mouse(slider.ax)
            except Exception:
                pass
            try:
                slider.disconnect_events()
            except Exception:
                pass
        self.alpha_sliders.clear()
        self.beta_sliders.clear()

    # -------------------------------------------------------------------
    # WINDOW 3: GEOMETRY CONTROLS  (blow-up/down)
    # -------------------------------------------------------------------

    def _build_control_figure(self):
        if self._ctrl_fig is not None:
            try:
                plt.close(self._ctrl_fig)
            except Exception:
                pass

        cone_lbls = self.geo_manager.cone_labels
        ray_lbls  = self.geo_manager.ray_labels
        n_cones   = len(cone_lbls)
        n_rays    = len(ray_lbls)

        fig_h = max(4.0, 0.35 * max(n_cones, n_rays) + 1.5)
        self._ctrl_fig = plt.figure("Geometry controls", figsize=(5.5, fig_h))
        self._ctrl_fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.08)

        cone_h    = 0.06 * n_cones
        ray_h     = 0.06 * n_rays
        radio_h   = max(cone_h, ray_h, 0.4)
        radio_top = 0.92

        ax_cone = self._ctrl_fig.add_axes([0.03, radio_top - radio_h, 0.45, radio_h])
        ax_cone.set_title("Blow-up cone", fontsize=9, pad=3)
        active_cone      = min(self._selected_cone, n_cones - 1)
        self._cone_radio = RadioButtons(ax_cone, cone_lbls, active=active_cone)
        self._cone_radio.on_clicked(self._on_cone_selected)

        ax_ray = self._ctrl_fig.add_axes([0.52, radio_top - radio_h, 0.45, radio_h])
        ax_ray.set_title("Blow-down ray", fontsize=9, pad=3)
        active_ray      = min(self._selected_ray, n_rays - 1)
        self._ray_radio = RadioButtons(ax_ray, ray_lbls, active=active_ray)
        self._ray_radio.on_clicked(self._on_ray_selected)

        btn_top = radio_top - radio_h - 0.04
        ax_up   = self._ctrl_fig.add_axes([0.10, btn_top - 0.12, 0.35, 0.10])
        ax_dn   = self._ctrl_fig.add_axes([0.55, btn_top - 0.12, 0.35, 0.10])
        self._blowup_btn   = Button(ax_up, "Blow up")
        self._blowdown_btn = Button(ax_dn, "Blow down")
        self._blowup_btn.on_clicked(self._on_blow_up)
        self._blowdown_btn.on_clicked(self._on_blow_down)

        self._ctrl_status = self._ctrl_fig.text(
            0.5, 0.02, "",
            ha="center", va="bottom", fontsize=9,
            transform=self._ctrl_fig.transFigure,
        )

        self._ctrl_fig.show()

    def _set_ctrl_status(self, msg: str):
        if self._ctrl_status is not None:
            self._ctrl_status.set_text(msg)
            self._ctrl_fig.canvas.draw_idle()

    # -------------------------------------------------------------------
    # WINDOW 4: INTERSECTION MATRIX
    # -------------------------------------------------------------------

    def _build_matrix_figure(self):
        if self._matrix_fig is not None:
            try:
                plt.close(self._matrix_fig)
            except Exception:
                pass

        M        = self.geo_manager.intersection_matrix
        n        = self.geo_manager.n_rays
        alphabet = string.ascii_lowercase

        if M is None:
            return

        fig_size         = max(3.0, 0.5 * n + 1.0)
        self._matrix_fig = plt.figure("Intersection matrix", figsize=(fig_size, fig_size))
        self._matrix_fig.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05)
        self._matrix_fig.suptitle("Intersection matrix", fontsize=10, fontstyle="italic")

        ax = self._matrix_fig.add_axes([0.05, 0.05, 0.90, 0.85])
        ax.axis("off")

        col_labels = [rf"$D_{{{alphabet[k]}}}$" for k in range(n)]
        row_labels = [rf"$D_{{{alphabet[k]}}}$" for k in range(n)]
        cell_text  = [[str(int(M[i, j])) for j in range(n)] for i in range(n)]

        cell_colors = []
        for i in range(n):
            row_colors = []
            for j in range(n):
                if i == j:
                    val = int(M[i, j])
                    if val == -1:
                        row_colors.append("#d4edda")
                    elif val < -1:
                        row_colors.append("#f8d7da")
                    else:
                        row_colors.append("#fff3cd")
                else:
                    row_colors.append("#ffffff")
            cell_colors.append(row_colors)

        table = ax.table(
            cellText=cell_text,
            rowLabels=row_labels,
            colLabels=col_labels,
            cellColours=cell_colors,
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.6)

        self._matrix_fig.show()

    # -------------------------------------------------------------------
    # GEOMETRY SETUP
    # -------------------------------------------------------------------

    def _apply_geometry(self):
        geo = self.geo_manager.current
        dim = geo.divisor_dimension

        self.alpha = np.ones(dim)
        self.beta  = np.ones(dim)

        max_k = dim - 1
        self.x_axis = AxisChoice(self.x_axis.vector, min(self.x_axis.k, max_k))
        self.y_axis = AxisChoice(self.y_axis.vector, min(self.y_axis.k, max_k))
        if self.x_axis == self.y_axis:
            self.y_axis = AxisChoice("beta", 0)

        self._selected_cone = 0
        self._selected_ray  = 0

        raw_ineq           = compute_nonredundant_inequalities(geo)
        self._ineq_strings = get_inequalities(raw_ineq)

        self._build_control_figure()
        self._build_controls_figure()
        self._build_sliders_figure()
        self._build_matrix_figure()

        self.plotter.im = None
        self.plotter.plot_fan(geo.rays)
        self.update_plot()

    # -------------------------------------------------------------------
    # SCHEDULE REBUILD  (axis change — rebuilds controls + sliders)
    # -------------------------------------------------------------------

    def _schedule_rebuild(self):
        if self._pending_rebuild:
            return
        self._pending_rebuild = True
        fig = self._controls_fig

        def _deferred(_):
            fig.canvas.mpl_disconnect(self._rebuild_cid)
            self._pending_rebuild = False
            self._remove_slider_widgets()
            self._build_controls_figure()
            self._build_sliders_figure()
            self.update_plot()

        self._rebuild_cid = fig.canvas.mpl_connect("draw_event", _deferred)
        fig.canvas.draw_idle()

    # -------------------------------------------------------------------
    # INEQUALITY COLORS
    # -------------------------------------------------------------------

    def _update_inequality_colors(self):
        if not self._ineq_texts:
            return

        alpha = self.alpha.copy()
        beta  = self.beta.copy()

        if self.x_axis.vector == "alpha":
            alpha[self.x_axis.k] = 0.5
        else:
            beta[self.x_axis.k]  = 0.5

        if self.y_axis.vector == "alpha":
            alpha[self.y_axis.k] = 0.5
        else:
            beta[self.y_axis.k]  = 0.5

        surface = self.geo_manager.current
        _, alpha_ints = is_kaehler(alpha, surface)
        _, beta_ints  = is_kaehler(beta,  surface)

        satisfied = [a > 0 and b > 0 for a, b in zip(alpha_ints, beta_ints)]

        for t, ok in zip(self._ineq_texts, satisfied):
            t.set_color("#2a7a2a" if ok else "#cc2222")

        if self._sliders_fig is not None:
            self._sliders_fig.canvas.draw()

    # -------------------------------------------------------------------
    # RADIO / BUTTON CALLBACKS
    # -------------------------------------------------------------------

    def _on_equation_changed(self, label):
        self.equation_name = label
        self.update_plot()
        self._update_inequality_colors()

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

    def _on_cone_selected(self, label):
        self._selected_cone = int(label.split(":")[0].split()[-1])

    def _on_ray_selected(self, label):
        self._selected_ray = int(label.split(":")[0].split()[-1])

    def _on_blow_up(self, _):
        try:
            self.geo_manager.blow_up(self._selected_cone)
        except Exception as e:
            self._set_ctrl_status(f"Blow-up failed: {e}")
            return
        self._apply_geometry()

    def _on_blow_down(self, _):
        try:
            self.geo_manager.blow_down(self._selected_ray)
        except Exception as e:
            self._set_ctrl_status(str(e))
            return
        self._apply_geometry()

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
    # PROGRAMMATIC AXIS SELECTION
    # -------------------------------------------------------------------

    def set_axes(self, x_axis: AxisChoice, y_axis: AxisChoice):
        if self._axes_conflict(x_axis, y_axis):
            raise ValueError("x_axis and y_axis must refer to different coefficients")
        self.x_axis = x_axis
        self.y_axis = y_axis
        self._remove_slider_widgets()
        self._build_controls_figure()
        self._build_sliders_figure()
        self.update_plot()

    # -------------------------------------------------------------------
    # COMPUTATION
    # -------------------------------------------------------------------

    def evaluate_slice(self):
        n          = len(self.axis_grid)
        X, Y       = np.meshgrid(self.axis_grid, self.axis_grid, indexing="ij")
        J          = np.zeros((n, n))
        mask       = np.zeros((n, n), dtype=bool)
        div_grid   = np.full((n, n), -1, dtype=int)
        base_alpha = self.alpha.copy()
        base_beta  = self.beta.copy()

        self.plotter.set_status("Computing...")
        self.plotter.fig.canvas.flush_events()

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

                val, valid, div_idx = evaluate_point(
                    self.equation_name,
                    (alpha, beta),
                    self.geo_manager.current,
                )
                J[a, b]    = val
                mask[a, b] = valid
                if div_idx is not None:
                    div_grid[a, b] = div_idx

        self.plotter.set_status("Done")
        return X, Y, J, mask, div_grid

    # -------------------------------------------------------------------
    # UPDATE LOOP
    # -------------------------------------------------------------------

    def update_plot(self):
        if self.geo_manager.current is None:
            return

        X, Y, J, mask, div_grid = self.evaluate_slice()

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

        self.plotter.plot_divisor_regions(X, Y, div_grid, mask)
        self.plotter.plot_zero_locus(X, Y, J, mask)
        self._update_inequality_colors()

    # -------------------------------------------------------------------
    # RUN
    # -------------------------------------------------------------------

    def run(self):
        group_name, variant_key = self._pick_surface()
        self.geo_manager.select_surface(group_name, variant_key)

        from plotting import Plotter
        self.plotter = Plotter()
        self.plotter.on_resolution_change = self._on_resolution_change

        self._apply_geometry()
        plt.show()

    def _on_resolution_change(self, new_res):
        self.resolution = new_res
        self.axis_grid  = np.linspace(np.finfo(float).eps, 1.0, new_res)
        self.update_plot()