import numpy as np
import sympy as sp
import inspect

from matplotlib.widgets import CheckButtons, Slider

# plt.subplot 객체의

x, a, b = sp.symbols('x a b')

def draw_function(expression, ax):
    """Plot an expression and its sampled maximum on [-1, 1]."""
    expression = sp.sympify(expression)
    if expression.free_symbols - {x}:
        raise ValueError("Substitute values for all symbols other than x before plotting.")
    x_values = np.linspace(-1, 1, 2000)
    function = sp.lambdify(x, expression, modules='numpy')
    y_values = np.broadcast_to(function(x_values), x_values.shape)
    line, = ax.plot(x_values, y_values)
    maximum = np.max(y_values)
    ax.axhline(maximum, color='red', linestyle='--',
               label=f'Maximum ≈ {maximum:.4g}')
    ax.set_xlabel('x')
    ax.grid(True)
    ax.legend()
    return line

def draw_basis_of_0val(F, ax):
    t_values = np.linspace(-1, 1, 1000)
    degree = len(F) - 1

    basis_funcs = []
    lines = []

    for i, F_i in enumerate(F):

        F_func = sp.lambdify(x, F_i, modules='numpy')

        basis_funcs.append(F_func)

        y_values = F_func(t_values)

        line, = ax.plot(t_values, y_values, label=rf'$F_{{{i}}}^{{{degree}}}(t)$')

        lines.append(line)
    
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel('t')
    ax.set_ylabel('Basis value')
    ax.grid(True)
    ax.legend()

    return t_values, basis_funcs, lines

def draw_basis_of_2vals(F, ax, a_init, b_init): #원하는 figure에 basis 그래프를 추가

    t_values = np.linspace(-1, 1, 1000)
    degree = len(F) - 1

    basis_funcs = []
    lines = []

    for i, F_i in enumerate(F):

        F_func = sp.lambdify((x,a,b), F_i, modules='numpy')

        basis_funcs.append(F_func)

        y_values = F_func(t_values, a_init, b_init)

        line, = ax.plot(t_values, y_values, label=rf'$F_{{{i}}}^{{{degree}}}(t)$')

        lines.append(line)
    
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel('t')
    ax.set_ylabel('Basis value')
    ax.grid(True)
    ax.legend()

    return t_values, basis_funcs, lines # 그래프를 그릴 t값, F_i(x,a)의 리스트, F_i의 Line2D 객체 리스트

def draw_basis_1val(F, ax, a_init):
    
    t_values = np.linspace(-1, 1, 1000)
    degree = len(F) - 1

    basis_funcs = []
    lines = []

    for i, F_i in enumerate(F):

        F_func = sp.lambdify((x,a), F_i, modules='numpy')

        basis_funcs.append(F_func)

        y_values = F_func(t_values, a_init)

        line, = ax.plot(t_values, y_values, label=rf'$F_{{{i}}}^{{{degree}}}(t)$')

        lines.append(line)
    
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel('t')
    ax.set_ylabel('Basis value')
    ax.grid(True)
    ax.legend()

    return t_values, basis_funcs, lines # 그래프를 그릴 t값, F_i(x,a)의 리스트, F_i의 Line2D 객체 리스트

def add_nodes(nodes, ax):
    return ax.scatter(nodes, np.zeros(len(nodes)))

def add_basis_checkboxes(basisArr, nodesArr, fig, ax):
    caller_locals = inspect.currentframe().f_back.f_locals

    labels = []

    for index, basis_set in enumerate(basisArr):
        label = next(
            (
                name for name, value in caller_locals.items()
                if value is basis_set
            ),
            f'basis {index + 1}'
        )
        labels.append(label)

    fig.subplots_adjust(right=0.78)

    checkbox_ax = fig.add_axes([0.82, 0.35, 0.15, 0.25])
    checkbox = CheckButtons(
        checkbox_ax,
        labels,
        [index == 0 for index in range(len(labels))]
    )

    # Matplotlib 이벤트는 위젯을 계속 보관하지 않으므로 축의 수명 동안 참조를 유지한다.
    checkbox_ax._basis_checkbox = checkbox

    updating = False

    def draw_selected_basis(selected_basis, selected_nodes):
        ax.clear()
        draw_basis_of_0val(selected_basis, ax)
        add_nodes(selected_nodes, ax)
        fig.canvas.draw_idle()

    def on_clicked(label):
        nonlocal updating

        if updating:
            return

        selected_index = labels.index(label)

        updating = True

        for index, status in enumerate(checkbox.get_status()):
            if index == selected_index:
                if not status:
                    checkbox.set_active(index)
            elif status:
                checkbox.set_active(index)

        updating = False

        draw_selected_basis(basisArr[selected_index], nodesArr[selected_index])

    checkbox.on_clicked(on_clicked)

    draw_selected_basis(basisArr[0], nodesArr[0])

    return checkbox

# 슬라이더 추가

def add_basis_slider(
    fig,
    ax,
    t_values,
    basis_funcs,
    lines,
    a_initial,
    scope
):
    
    fig.subplots_adjust(bottom=0.25)

    ax_slider = fig.add_axes(
        [0.2, 0.1, 0.65, 0.04]
    )

    a_slider = Slider(
        ax=ax_slider,
        label='a',
        valmin=scope[0],
        valmax=scope[1],
        valinit=scope[2],
        valstep=scope[3]
    )

    def update(value):
        current_a = value

        for F_func, line in zip(
            basis_funcs,
            lines
        ):
            new_y = F_func(
                t_values,
                current_a
            )

            line.set_ydata(new_y)

        ax.relim()
        ax.autoscale_view()

        ax.set_title(
            rf'Basis functions for $a={current_a:.2f}$'
        )

        fig.canvas.draw_idle()

    def on_key_press(event):
        step = a_slider.valstep

        if event.key == 'right':
            new_value = a_slider.val + step

        elif event.key == 'left':
            new_value = a_slider.val - step

        else:
            return

        new_value = max(
            a_slider.valmin,
            min(
                a_slider.valmax,
                new_value
            )
        )

        a_slider.set_val(new_value)

    a_slider.on_changed(update)

    fig.canvas.mpl_connect(
        'key_press_event',
        on_key_press
    )

    return a_slider

def add_basis_sliders(
    fig,
    ax,
    t_values,
    basis_funcs,
    lines,
    a_initial,
    b_initial,
    scope_a,
    scope_b
):
    # 슬라이더 두 개를 위한 공간
    fig.subplots_adjust(bottom=0.32)

    # a 슬라이더 축
    ax_slider_a = fig.add_axes(
        [0.2, 0.15, 0.65, 0.04]
    )

    # b 슬라이더 축
    ax_slider_b = fig.add_axes(
        [0.2, 0.08, 0.65, 0.04]
    )

    a_slider = Slider(
        ax=ax_slider_a,
        label='a',
        valmin=scope_a[0],
        valmax=scope_a[1],
        valinit=scope_a[2],
        valstep=scope_a[3]
    )

    b_slider = Slider(
        ax=ax_slider_b,
        label='b',
        valmin=scope_b[0],
        valmax=scope_b[1],
        valinit=scope_b[2],
        valstep=scope_b[3]
    )

    # 유효한 마지막 슬라이더 값 저장
    previous_values = {
        'a': a_initial,
        'b': b_initial
    }

    # 콜백 함수가 서로 다시 호출되는 것을 방지
    updating = False

    def update(value):
        nonlocal updating

        if updating:
            return

        current_a = a_slider.val
        current_b = b_slider.val

        # 반드시 0 < b < a < 1이어야 함
        if current_b >= current_a:
            updating = True

            # 조건을 위반하면 직전의 유효한 값으로 복원
            a_slider.set_val(previous_values['a'])
            b_slider.set_val(previous_values['b'])

            updating = False
            return

        previous_values['a'] = current_a
        previous_values['b'] = current_b

        # 모든 기저함수의 y값 갱신
        for F_func, line in zip(
            basis_funcs,
            lines
        ):
            new_y = F_func(
                t_values,
                current_a,
                current_b
            )

            line.set_ydata(new_y)

        # 새 데이터에 맞게 y축 범위 조절
        ax.relim()
        ax.autoscale_view()

        ax.set_title(
            rf'$a={current_a:.2f},\quad '
            rf'b={current_b:.2f}$'
        )

        fig.canvas.draw_idle()

    a_slider.on_changed(update)
    b_slider.on_changed(update)

    def on_key_press(event):
        step = scope_a[3]

        current_a = a_slider.val
        current_b = b_slider.val

        # 좌우 방향키: a 조절
        if event.key == 'right':
            new_a = min(
                a_slider.valmax,
                current_a + step
            )
            a_slider.set_val(new_a)

        elif event.key == 'left':
            new_a = max(
                current_b + step,
                current_a - step
            )
            a_slider.set_val(new_a)

        # 위아래 방향키: b 조절
        elif event.key == 'up':
            new_b = min(
                current_a - step,
                current_b + step
            )
            b_slider.set_val(new_b)

        elif event.key == 'down':
            new_b = max(
                b_slider.valmin,
                current_b - step
            )
            b_slider.set_val(new_b)

    fig.canvas.mpl_connect(
        'key_press_event',
        on_key_press
    )
    
    return a_slider, b_slider
    


def draw_basis_selector(sources, n_initial=4, n_min=2, n_max=12, initial_basis=None):
    """기저 선택과 차수 선택을 한 figure로 제공한다.

    sources는 {이름: 노드 생성 함수 또는 custom 노드} 딕셔너리다.
    generate_nodes 모듈에 등록된 함수만 가변 n을 지원한다.
    custom은 노드 배열(길이와 같은 n에서만 표시) 또는 {n: 노드 배열}로 받는다.
    custom에 해당 n이 없으면 그래프를 비우되 선택 상태는 유지한다.
    선택은 항상 하나이며 (fig, ax, radio, slider)를 반환한다.
    선택한 조합만 계산하고 최근 32개 조합의 그래프 데이터를 캐시한다.
    슬라이더 클릭/좌우 방향키를 지원하며 plt.show()는 호출자가 실행한다.
    """
    from functools import lru_cache
    from matplotlib.widgets import RadioButtons
    import matplotlib.pyplot as plt
    import generate_basis as basis_generator
    import generate_nodes as node_generators

    if not isinstance(sources, dict) or not sources:
        raise ValueError("Provide a non-empty dictionary of named node sources.")
    if any(not isinstance(name, str) or not name.strip() for name in sources):
        raise ValueError("Basis names must be non-empty strings.")
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer))
           for v in (n_initial, n_min, n_max)):
        raise ValueError("Degree values must be integers.")
    if not 1 <= n_min <= n_initial <= n_max or n_min == n_max:
        raise ValueError("Require 1 <= n_min <= n_initial <= n_max and n_min < n_max.")
    labels = list(sources)
    selected = labels[0] if initial_basis is None else initial_basis
    if selected not in sources:
        raise ValueError("initial_basis must be a name in sources.")

    def validate_nodes(raw, degree):
        nodes = np.asarray(raw)
        if np.iscomplexobj(nodes):
            if np.any(np.abs(nodes.imag) > 1e-10):
                raise ValueError("Nodes must be real.")
            nodes = nodes.real
        nodes = np.asarray(nodes, dtype=float)
        if (nodes.ndim != 1 or len(nodes) != degree or not np.all(np.isfinite(nodes))
                or len(np.unique(nodes)) != degree):
            raise ValueError("Expected n distinct finite nodes.")
        return nodes.copy()

    registered = [value for value in vars(node_generators).values() if inspect.isfunction(value)
                  and value.__module__ == node_generators.__name__]
    normalized = {}
    for name, source in sources.items():
        if callable(source):
            if not any(source is function for function in registered):
                raise ValueError(f"{name}: use a function defined in generate_nodes.py.")
            normalized[name] = source
        else:
            mapping = source if isinstance(source, dict) else {len(source): source}
            custom = {}
            for degree, raw in mapping.items():
                if isinstance(degree, (bool, np.bool_)) or not isinstance(degree, (int, np.integer)) or degree < 1:
                    raise ValueError("Custom degree keys must be positive integers.")
                custom[int(degree)] = validate_nodes(raw, degree)
            normalized[name] = custom
    t = np.linspace(-1, 1, 1000)

    @lru_cache(maxsize=32)
    def prepare(name, degree):
        source = normalized[name]
        if not callable(source) and degree not in source:
            return None
        nodes = validate_nodes(source(degree) if callable(source) else source[degree], degree)
        F = basis_generator.get_New_Basis(nodes, basis_name=name)
        values = np.array([np.broadcast_to(sp.lambdify(x, f, 'numpy')(t), t.shape)
                           for f in F], dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("Non-finite basis values; check normalization weights.")
        return nodes, values

    initial = prepare(selected, int(n_initial))
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(bottom=0.25, right=0.76)
    radio_ax = fig.add_axes([0.79, 0.30, 0.19, 0.55])
    radio = RadioButtons(radio_ax, labels, active=labels.index(selected))
    slider_ax = fig.add_axes([0.2, 0.10, 0.5, 0.04])
    slider = Slider(slider_ax, 'n', n_min, n_max, valinit=n_initial,
                    valstep=1, valfmt='%d', dragging=False)
    status = fig.text(0.1, 0.02, '', color='firebrick', fontsize=9)
    current_n = int(n_initial)
    updating = False

    def render(name, degree, data):
        ax.clear()
        ax.set(xlabel='t', ylabel='Basis value', xlim=(-1, 1), title=f'{name}: n={degree}')
        if data is None:
            # custom에 정의되지 않은 차수: 이전 선과 노드를 모두 제거한다.
            return
        nodes, values = data
        for index, y in enumerate(values):
            ax.plot(t, y, label=rf'$F_{{{index}}}^{{{degree}}}(t)$')
        add_nodes(nodes, ax)
        ax.axhline(0, linewidth=0.8)
        ax.grid(True)
        ax.legend(loc='upper right', ncol=2 if len(values) > 8 else 1)

    def update(name, degree):
        nonlocal selected, current_n, updating
        if updating or (name == selected and degree == current_n):
            return
        try:
            data = prepare(name, degree)
        except Exception as error:
            status.set_text(f'{name}, n={degree}: {error}')
            updating = True
            try:
                radio.set_active(labels.index(selected))
                slider.set_val(current_n)
            finally:
                updating = False
        else:
            selected, current_n = name, degree
            render(name, degree, data)
            status.set_text('')
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key in ('left', 'right'):
            step = 1 if event.key == 'right' else -1
            slider.set_val(max(n_min, min(n_max, current_n + step)))

    radio.on_clicked(lambda name: update(name, current_n))
    slider.on_changed(lambda value: update(selected, int(value)))
    fig.canvas.mpl_connect('key_press_event', on_key)
    # figure가 살아 있는 동안 모든 위젯과 콜백을 보관한다.
    fig._basis_selector_controls = (radio, slider, on_key)
    render(selected, current_n, initial)
    return fig, ax, radio, slider


def draw_adjustable_generalized_lebesgue_function(basis, a_initial=0.5, a_min=0.01, a_max=0.99, step=0.01):
    """매개변수 a가 포함된 기저의 일반화 르베그 함수를 슬라이더로 살펴본다."""
    import matplotlib.pyplot as plt
    from generalized_Lebegue_constant import generalized_Lebegue_function

    if not 0 < a_min <= a_initial <= a_max < 1 or step <= 0 or a_min == a_max:
        raise ValueError('Require 0 < a_min <= a_initial <= a_max < 1 and step > 0.')
    expression = generalized_Lebegue_function(basis)
    if expression.free_symbols - {x, a}:
        raise ValueError('The basis may contain only x and a.')
    function = sp.lambdify((x, a), expression, modules='numpy')
    x_values = np.linspace(-1, 1, 1000)
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(bottom=0.22)
    slider = Slider(fig.add_axes([0.2, 0.08, 0.6, 0.04]), 'a', a_min, a_max,
                    valinit=a_initial, valstep=step)
    line, = ax.plot([], [])
    maximum_line = ax.axhline(0, color='red', linestyle='--')
    ax.set(xlabel='x', ylabel=r'$\Lambda_F(x)$', xlim=(-1, 1))
    ax.grid(True)

    n = len(basis)
    def update(value):
        values = np.asarray(np.broadcast_to(function(x_values, value), x_values.shape), dtype=float)
        line.set_data(x_values, values)
        maximum = np.max(values)
        maximum_line.set_ydata([maximum, maximum])
        maximum_line.set_label(f'Maximum ≈ {maximum:.4g}')
        ax.set_title(f'n={n-1}, Generalized Lebesgue function: a={value:.2f}')
        ax.relim()
        ax.autoscale_view(scalex=False)
        ax.legend(loc="upper right")
        fig.canvas.draw_idle()

    slider.on_changed(update)
    fig._adjustable_generalized_lebesgue_function_controls = (slider, update)
    update(a_initial)
    return fig, ax, slider


def draw_generalized_lebesgue_function_selector(sources, n_initial=4, n_min=2, n_max=12, initial_basis=None):
    """체크박스로 여러 기저의 일반화 르베그 함수를 비교하고 좌우 키로 n을 조절한다."""
    from functools import lru_cache
    import matplotlib.pyplot as plt
    import generate_basis as basis_generator
    from generalized_Lebegue_constant import generalized_Lebegue_function

    if not isinstance(sources, dict) or not sources or any(
        not isinstance(name, str) or not name.strip() for name in sources
    ):
        raise ValueError('Provide a non-empty dictionary of named node sources.')
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer))
           for v in (n_initial, n_min, n_max)) or not 1 <= n_min <= n_initial <= n_max or n_min == n_max:
        raise ValueError('Require integer values 1 <= n_min <= n_initial <= n_max and n_min < n_max.')
    labels = list(sources)
    initial = labels[0] if initial_basis is None else initial_basis
    if initial not in sources:
        raise ValueError('initial_basis must be a name in sources.')
    x_values = np.linspace(-1, 1, 1000)

    @lru_cache(maxsize=32)
    def prepare(name, degree):
        source = sources[name]
        if callable(source):
            raw = source(degree)
        else:
            mapping = source if isinstance(source, dict) else {len(source): source}
            if degree not in mapping:
                return None
            raw = mapping[degree]
        nodes = np.asarray(raw)
        if np.iscomplexobj(nodes) and np.any(np.abs(nodes.imag) > 1e-10):
            raise ValueError('Nodes must be real.')
        nodes = np.asarray(nodes.real, dtype=float)
        if nodes.ndim != 1 or len(nodes) != degree or not np.all(np.isfinite(nodes)) or len(np.unique(nodes)) != degree:
            raise ValueError('Expected n distinct finite nodes.')
        expression = generalized_Lebegue_function(basis_generator.get_New_Basis(nodes, basis_name=name))
        values = np.asarray(np.broadcast_to(sp.lambdify(x, expression, 'numpy')(x_values), x_values.shape), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError('Non-finite values; check normalization weights.')
        return values

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(bottom=0.25, right=0.76)
    if 'ALL' in labels:
        raise ValueError("'ALL' is reserved for the select-all checkbox.")
    checks = CheckButtons(fig.add_axes([0.79, 0.52, 0.19, 0.38]), ['ALL', *labels],
                          [len(labels) == 1, *[name == initial for name in labels]])
    max_ax = fig.add_axes([0.79, 0.16, 0.19, 0.30])
    max_ax.set_xticks([])
    max_ax.set_yticks([])
    max_check = CheckButtons(fig.add_axes([0.80, 0.40, 0.16, 0.05]), ['Max'], [False])
    max_text = max_ax.text(0.08, 0.77, '', transform=max_ax.transAxes,
                           va='top', fontsize=9)
    slider = Slider(fig.add_axes([0.2, 0.10, 0.5, 0.04]), 'n', n_min, n_max,
                    valinit=n_initial, valstep=1, valfmt='%d', dragging=False)
    status = fig.text(0.1, 0.02, '', color='firebrick', fontsize=9)

    def render(_=None):
        degree = int(slider.val)
        ax.clear()
        ax.set(xlabel='x', ylabel=r'$\Lambda_F(x)$', xlim=(-1, 1),
               title=f'Generalized Lebesgue functions: n={degree}')
        errors = []
        maxima = []
        for index, (name, active) in enumerate(zip(labels, checks.get_status()[1:])):
            if not active:
                continue
            try:
                values = prepare(name, degree)
                if values is None:
                    errors.append(f'{name}: no nodes for n={degree}')
                    continue
                ax.plot(x_values, values, color=f'C{index % 10}', label=name)
                maximum = np.max(values)
                maxima.append(f'{name}: {maximum:.4g}')
                if max_check.get_status()[0]:
                    ax.axhline(maximum, color=f'C{index % 10}', linestyle='--')
            except Exception as error:
                errors.append(f'{name}: {error}')
        ax.grid(True)
        if ax.lines:
            ax.legend()
        status.set_text('; '.join(errors))
        max_text.set_text('\n'.join(maxima))
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key in ('left', 'right'):
            step = 1 if event.key == 'right' else -1
            slider.set_val(max(n_min, min(n_max, int(slider.val) + step)))

    updating_checks = False

    def on_check(label):
        nonlocal updating_checks
        if updating_checks:
            return
        updating_checks = True
        try:
            states = checks.get_status()
            if label == 'ALL':
                target = not all(states[1:])
                for index, active in enumerate(states[1:], start=1):
                    if active != target:
                        checks.set_active(index)
            all_selected = all(checks.get_status()[1:])
            if checks.get_status()[0] != all_selected:
                checks.set_active(0)
        finally:
            updating_checks = False
        render()

    checks.on_clicked(on_check)
    max_check.on_clicked(render)
    slider.on_changed(render)
    fig.canvas.mpl_connect('key_press_event', on_key)
    fig._generalized_lebesgue_function_selector_controls = (checks, slider, on_key, max_check)
    render()
    return fig, ax, checks, slider


def draw_generalized_lebesgue_constant(sources, n_initial=4, n_min=2, n_max=12, initial_basis=None):
    """선택한 기저의 상수를 n_min부터 선택한 n까지 비교한다.

    sources는 {이름: 노드 생성 함수 또는 사용자 노드 배열/차수별 배열}이다.
    체크박스와 좌우 방향키 또는 슬라이더로 선택한다. 상수는 4001개 x 표본의
    근삿값이며, 반환값은 (figure, axes, 체크박스, 슬라이더)이다.
    """
    from functools import lru_cache
    import matplotlib.pyplot as plt
    from generalized_Lebegue_constant import generalized_constant_for_nodes

    if not isinstance(sources, dict) or not sources or any(
            not isinstance(name, str) or not name.strip() for name in sources):
        raise ValueError('Provide a non-empty dictionary of named node sources.')
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer))
           for v in (n_initial, n_min, n_max)) or not 1 <= n_min <= n_initial <= n_max or n_min == n_max:
        raise ValueError('Require integer values 1 <= n_min <= n_initial <= n_max and n_min < n_max.')
    labels = list(sources)
    if 'ALL' in labels:
        raise ValueError("'ALL' is reserved for the select-all checkbox.")
    initial = labels[0] if initial_basis is None else initial_basis
    if initial not in sources:
        raise ValueError('initial_basis must be a name in sources.')

    @lru_cache(maxsize=256)
    def prepare(name, degree):
        source = sources[name]
        if callable(source):
            raw = source(degree)
        else:
            mapping = source if isinstance(source, dict) else {len(source): source}
            if degree not in mapping:
                return None
            raw = mapping[degree]
        nodes = np.asarray(raw)
        if nodes.ndim != 1 or len(nodes) != degree:
            raise ValueError('Expected n nodes.')
        return generalized_constant_for_nodes(nodes, basis_name=name)

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(bottom=0.25, right=0.76)
    checks = CheckButtons(fig.add_axes([0.79, 0.52, 0.19, 0.38]), ['ALL', *labels],
                          [len(labels) == 1, *[name == initial for name in labels]])
    info_ax = fig.add_axes([0.79, 0.16, 0.19, 0.30])
    info_ax.set_xticks([])
    info_ax.set_yticks([])
    info = info_ax.text(0.08, 0.92, '', transform=info_ax.transAxes, va='top', fontsize=9)
    slider = Slider(fig.add_axes([0.2, 0.10, 0.5, 0.04]), 'n', n_min, n_max,
                    valinit=n_initial, valstep=1, valfmt='%d', dragging=False)
    status = fig.text(0.1, 0.02, '', color='firebrick', fontsize=9)

    def render(_=None):
        degree = int(slider.val)
        ax.clear()
        ax.set(xlabel='n (number of nodes)', ylabel='Generalized Lebesgue constant (approx.)',
               title=f'Generalized Lebesgue constants: n={degree}')
        ax.set_xticks(range(n_min, degree + 1))
        ax.set_xlim(n_min - 0.5, degree + 0.5)
        errors, current = [], []
        for index, (name, active) in enumerate(zip(labels, checks.get_status()[1:])):
            if not active:
                continue
            degrees, values = [], []
            for n in range(n_min, degree + 1):
                try:
                    value = prepare(name, n)
                except Exception as error:
                    errors.append(f'{name}, n={n}: {error}')
                    value = None
                degrees.append(n)
                values.append(np.nan if value is None else value)
            if np.any(np.isfinite(values)):
                ax.plot(degrees, values, marker='o', color=f'C{index % 10}', label=name)
            current.append(f'{name}: {values[-1]:.4g}' if np.isfinite(values[-1]) else f'{name}: unavailable')
        ax.grid(True)
        if ax.lines:
            ax.legend()
        info.set_text(f'n = {degree}\n' + '\n'.join(current))
        status.set_text('; '.join(errors))
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key in ('left', 'right'):
            step = 1 if event.key == 'right' else -1
            slider.set_val(max(n_min, min(n_max, int(slider.val) + step)))

    updating_checks = False

    def on_check(label):
        nonlocal updating_checks
        if updating_checks:
            return
        updating_checks = True
        try:
            states = checks.get_status()
            if label == 'ALL':
                target = not all(states[1:])
                for index, active in enumerate(states[1:], start=1):
                    if active != target:
                        checks.set_active(index)
            all_selected = all(checks.get_status()[1:])
            if checks.get_status()[0] != all_selected:
                checks.set_active(0)
        finally:
            updating_checks = False
        render()

    checks.on_clicked(on_check)
    slider.on_changed(render)
    fig.canvas.mpl_connect('key_press_event', on_key)
    fig._generalized_lebesgue_constant_controls = (checks, slider, on_key)
    render()
    return fig, ax, checks, slider


def draw_basis_by_degree(node_function, n_initial=4, n_min=2, n_max=12):
    """단일 노드 생성 함수용 기존 인터페이스. 통합 선택기를 재사용한다."""
    name = getattr(node_function, '__name__', 'basis')
    fig, ax, _, slider = draw_basis_selector(
        {name: node_function}, n_initial, n_min, n_max
    )
    return fig, ax, slider
