import numpy as np

f = np.sin
df = np.cos
x0 = 1.0
exact = df(x0)

print(f"exact value         = {exact:.10f}")

for h in [0.1, 0.01, 0.001, 0.0001]:
    d = (f(x0 + h) - f(x0)) / h
    print(f"h ={h:6.0e}: approximation = {d:.10f}, error = {abs(d - exact):.3e}")
print()

for h in [0.1, 0.01, 0.001, 0.0001]:
    d = (f(x0 + h) - f(x0 - h)) / (2*h)
    print(f"h ={h:6.0e}: approximation = {d:.10f}, error = {abs(d - exact):.3e}")
print()


import matplotlib.pyplot as plt

hs = np.logspace(-15, -1, 400)
err_f = np.abs((f(x0 + hs) - f(x0)) / hs - exact)
err_c = np.abs((f(x0 + hs) - f(x0 - hs)) / (2 * hs) - exact)

fig, ax = plt.subplots(figsize=(8, 5))
ax.loglog(hs, err_f, lw=2, label="forward difference")
ax.loglog(hs, err_c, lw=2, label="central difference")
ax.set_xlabel("h")
ax.set_ylabel("error")
ax.grid(alpha=0.3, which="both")
ax.legend()
plt.show()

print("오차가 가장 작은 h =", hs[np.argmin(err_c)])

# ------------------------------------------------------------
# mpmath : 정밀도(자릿수)를 늘리면 round-off error 가 줄어든다
# ------------------------------------------------------------
# 전체 오차 ≈ 절단오차 (h^2/6)|f'''| + 반올림오차 (eps/h)
#   - 절단오차는 정밀도와 무관하다 (수식 자체의 오차)
#   - 반올림오차의 eps 는 자릿수(dps)가 늘수록 10^(-dps) 로 작아진다
# 그래서 dps 를 늘리면 V 자 그래프의 오른쪽(절단오차)은 그대로이고,
# 왼쪽(반올림오차)만 아래로 내려간다.
from mpmath import mp, mpf

ks = range(1, 40)                     # h = 10^-1, 10^-2, ..., 10^-39

fig, ax = plt.subplots(figsize=(8, 5))
for dps in [15, 30, 60]:              # 15자리 ≈ float64
    mp.dps = dps
    x = mpf(1)
    hs_mp, errs = [], []
    for k in ks:
        h = mpf(10) ** (-k)
        d = (mp.sin(x + h) - mp.sin(x - h)) / (2 * h)
        hs_mp.append(h)
        errs.append(abs(d - mp.cos(x)))
    best = min(range(len(errs)), key=lambda i: errs[i])
    print(f"dps = {dps:2d}: 가장 좋은 h = 1e-{ks[best]}, 최소 오차 = {mp.nstr(errs[best], 3)}")
    ax.loglog([float(h) for h in hs_mp], [float(e) for e in errs], "o-", label=f"mpmath dps = {dps}")

ax.set_xlabel("h")
ax.set_ylabel("error")
ax.set_title("central difference : round-off error vs precision")
ax.grid(alpha=0.3, which="both")
ax.legend()
plt.show()
