import numpy as np
import sympy as sp

from matplotlib.widgets import CheckButtons

x = sp.symbols('x')

# basis 배열과 조절점 배열, 그래프 
def draw_interactive_curves(basis_sets, P, ax, labels=None):
    t_values = np.linspace(-1, 1, 1000)

    px_values = list(P[0])
    py_values = list(P[1])

    if labels is None:
        labels = [f'curve {i + 1}' for i in range(len(basis_sets))]

    basis_functions = [
        [sp.lambdify(x, basis, modules='numpy') for basis in basis_set]
        for basis_set in basis_sets
    ]

    def evaluate_curve(basis_func_set):
        curve_x_values = np.zeros_like(t_values, dtype=float)
        curve_y_values = np.zeros_like(t_values, dtype=float)

        for i, basis_func in enumerate(basis_func_set):
            basis_values = basis_func(t_values)
            curve_x_values += px_values[i] * basis_values
            curve_y_values += py_values[i] * basis_values

        return curve_x_values, curve_y_values

    curve_lines = []

    for basis_func_set, label in zip(basis_functions, labels):
        curve_x_values, curve_y_values = evaluate_curve(basis_func_set)
        line, = ax.plot(
            curve_x_values,
            curve_y_values,
            linewidth=2,
            label=label
        )
        curve_lines.append((line, basis_func_set))

    control_polygon_line, = ax.plot(
        px_values,
        py_values,
        'k--',
        linewidth=1,
        label='control polygon'
    )

    control_point_markers = []
    control_point_labels = []

    for i, (point_x, point_y) in enumerate(zip(px_values, py_values)):
        marker, = ax.plot(point_x, point_y, 'ro', markersize=7)
        label = ax.text(point_x, point_y, f'P{i}')
        control_point_markers.append(marker)
        control_point_labels.append(label)

    def update_control_points():
        for i, marker in enumerate(control_point_markers):
            marker.set_data([px_values[i]], [py_values[i]])
            control_point_labels[i].set_position((px_values[i], py_values[i]))

    def find_control_point_index(event):
        if event.xdata is None or event.ydata is None:
            return None

        threshold_pixels = 12

        for i, (point_x, point_y) in enumerate(zip(px_values, py_values)):
            display_x, display_y = ax.transData.transform((point_x, point_y))
            distance = ((event.x - display_x) ** 2 + (event.y - display_y) ** 2) ** 0.5

            if distance <= threshold_pixels:
                return i

        return None

    def redraw_curve():
        for curve_line, basis_func_set in curve_lines:
            curve_x_values, curve_y_values = evaluate_curve(basis_func_set)
            curve_line.set_data(curve_x_values, curve_y_values)

        control_polygon_line.set_data(px_values, py_values)
        update_control_points()
        ax.relim()
        ax.autoscale_view()
        ax.axis('equal')
        fig.canvas.draw_idle()

    fig = ax.figure
    dragged_point_index = None

    def on_press(event):
        nonlocal dragged_point_index

        if event.inaxes != ax:
            return

        dragged_point_index = find_control_point_index(event)

    def on_motion(event):
        if dragged_point_index is None or event.inaxes != ax:
            return

        if event.xdata is None or event.ydata is None:
            return

        px_values[dragged_point_index] = event.xdata
        py_values[dragged_point_index] = event.ydata
        redraw_curve()

    def on_release(_event):
        nonlocal dragged_point_index
        dragged_point_index = None

    fig.canvas.mpl_connect('button_press_event', on_press)
    fig.canvas.mpl_connect('motion_notify_event', on_motion)
    fig.canvas.mpl_connect('button_release_event', on_release)

    ax.set_xlabel('x(t)')
    ax.set_ylabel('y(t)')
    ax.axis('equal')
    ax.grid(True)
    ax.legend()

    return [line for line, _basis_func_set in curve_lines]


def add_curve_checkboxes(labels, lines, fig, ax):
    fig.subplots_adjust(right=0.78)

    checkbox_ax = fig.add_axes([0.82, 0.35, 0.15, 0.25])
    checkbox = CheckButtons(
        checkbox_ax,
        labels,
        [line.get_visible() for line in lines]
    )

    checkbox_ax._curve_checkbox = checkbox

    def refresh_legend():
        visible_lines = [
            line for line in lines
            if line.get_visible()
        ]

        if visible_lines:
            ax.legend(handles=visible_lines)
        else:
            legend = ax.get_legend()

            if legend is not None:
                legend.remove()

    def on_clicked(label):
        selected_index = labels.index(label)
        selected_line = lines[selected_index]

        selected_line.set_visible(
            not selected_line.get_visible()
        )

        refresh_legend()
        fig.canvas.draw_idle()

    checkbox.on_clicked(on_clicked)

    return checkbox
