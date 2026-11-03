# 임의 개수의 조절 변수

`draw_adjustable_basis`에 기호가 포함된 노드와 각 기호의 `(최솟값, 최댓값, 초기값, 간격)`을 전달하면 변수마다 슬라이더가 생성됩니다.

```python
import sympy as sp
import matplotlib.pyplot as plt
import draw_basis as draw

a, b, c = sp.symbols('a b c')
fig, ax, sliders = draw.draw_adjustable_basis(
    [-1, a, b, c, 1],
    {a: (-0.9, -0.1, -0.5, 0.01),
     b: (-0.4, 0.4, 0.0, 0.01),
     c: (0.1, 0.9, 0.5, 0.01)},
)
plt.show()
```

노드가 겹치거나 기저가 정의되지 않는 값으로 이동하면 마지막 유효 값으로 복원합니다. 변수 개수가 늘면 기호 계산에 시간이 오래 걸릴 수 있습니다. 일반화 르베그 함수와 상수의 기존 그래프는 별도의 단일 변수 `a` 예제입니다.
