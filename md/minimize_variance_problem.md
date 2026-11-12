<h1 style="background-color: rgb(134, 227, 182); color: rgb(0, 0, 0)"> Critical Point Problem </h1>

# Definitions

Let $\{\tau_1, \tau_2, \dots, \tau_n\}$ be arbitrary distinct ordered nodes with $\tau_i \in [-1, 1]$.

### Lagrange polynomial $L_i^n(t)$ associated with $\tau_i$

$$
L_i^n(t) = \prod_{\substack{j=1\\j\ne i}}^n \frac{t-\tau_j}{\tau_i-\tau_j}, \qquad i=1,\dots,n.
$$

For convenience, define $D_i^n = \prod_{\substack{j=1\\j\ne i}}^n (\tau_i-\tau_j)$.



### Quadrature weight $\omega_i^n$ associated with $\tau_i$

$$
\omega_i^n = \int_{-1}^{1} L_i^n(t)\,dt, \qquad i=1,\dots,n.
$$

### Normalized Lagrange polynomial $\tilde L_i^n(t)$

Throughout, we assume $\omega_i^n\ne0$ for all $i=1,\dots,n$, so that $\tilde L_i^n$ is well defined. (For arbitrary nodes the weights may be negative, but they must not vanish.)

$$
\tilde L_i^n(t) = \frac{L_i^n(t)}{\omega_i^n}, \qquad i=1,\dots,n.
$$

### Cumulative Lagrange polynomial $G_i^n(t)$

$$
G_i^n(t)=\int_{-1}^{t}\tilde L_i^n(s)\,ds-\frac12, \qquad i=1,\dots,n.
$$

For the endpoint indices, define $G_0^n(t)=\frac12$ and $G_{n+1}^n(t)=-\frac12$.

### Basis polynomial $F_i^n(t)$

$$
F_i^n(t)=G_i^n(t)-G_{i+1}^n(t), \qquad i=0,1,\dots,n.
$$

For $i=1,\dots,n-1$, the constants $-\frac12$ cancel, so

$$
\begin{aligned}
F_i^n(t)&=\int_{-1}^{t}\frac{L_i^n(s)}{\omega_i^n}\,ds-\int_{-1}^{t}\frac{L_{i+1}^n(s)}{\omega_{i+1}^n}\,ds\\
&=\int_{-1}^{t}\left(\frac{L_i^n(s)}{\omega_i^n}-\frac{L_{i+1}^n(s)}{\omega_{i+1}^n}\right)ds,
\qquad i=1,\dots,n-1.
\end{aligned}
$$

For the endpoint indices,

$$
F_0^n(t)=1-\int_{-1}^{t}\tilde L_1^n(s)\,ds,\qquad
F_n^n(t)=\int_{-1}^{t}\tilde L_n^n(s)\,ds.
$$

<br><br><br>

# Critical points of $F_i^n(t)$ with the ordered nodes $\{\tau_1,\tau_2,\dots,\tau_n\}$

## 1) $i=1,\dots,n-1$

$$
\begin{aligned}
\frac{d}{dt}F_i^n(t)=0
&\iff \frac{L_i^n(t)}{\omega_i^n}=\frac{L_{i+1}^n(t)}{\omega_{i+1}^n}\\
&\iff \prod_{\substack{j=1\\j\notin\{i,i+1\}}}^{n}(t-\tau_j)\left\{\alpha(t-\tau_{i+1})-\beta(t-\tau_i)\right\}=0,
\end{aligned}
$$

where $\alpha=\frac{\omega_{i+1}^n}{D_i^n}$ and $\beta=\frac{\omega_i^n}{D_{i+1}^n}$.

If $\alpha\ne\beta$, the critical points are $t=\tau_k$ for $k\ne i,i+1$, or

$$
t=t^*:=\frac{\alpha\tau_{i+1}-\beta\tau_i}{\alpha-\beta}.
$$

## 2) $i=0,n$

$$
\frac{d}{dt}F_0^n(t)=-\tilde L_1^n(t)\quad\Rightarrow\quad t=\tau_2,\tau_3,\dots,\tau_n,
$$

$$
\frac{d}{dt}F_n^n(t)=\tilde L_n^n(t)\quad\Rightarrow\quad t=\tau_1,\tau_2,\dots,\tau_{n-1}.
$$

There does not exist $t^*$ at the endpoint basis indices.

<br><br><br>

# Special cases of $t^*$

## 1) $\alpha=\beta$

In this case,

$$
(\alpha-\beta)t=\alpha\tau_{i+1}-\beta\tau_i=\alpha(\tau_{i+1}-\tau_i)\ne0,
$$

so the linear equation determining the additional critical point has no real solution. Thus $F$ only has $n-2$ critical points.

## 2) $t^*=\tau_k$, $k\ne i$, $k\ne i+1$

In this case, $\frac{d}{dt}F$ has a double root at $\tau_k$.

<br><br><br>

# Sign of $\frac{d^2}{dt^2}F_i^n(\tau_k)$, $k\ne i$, $k\ne i+1$

Assume $\omega_j^n>0$ for all $j=1,\dots,n$ throughout this section.

Consider

$$
\frac{d}{dt}\tilde L_i^n(t)=\frac{1}{\omega_i^nD_i^n}\frac{d}{dt}\left(\prod_{\substack{j=1\\j\ne i}}^n(t-\tau_j)\right).
$$

$$
\frac{d}{dt}\tilde L_i^n(\tau_k)=\frac{1}{\omega_i^nD_i^n}\prod_{\substack{j=1\\j\ne i,k}}^n(\tau_k-\tau_j).
$$

Since the nodes are ordered, the sign of the product is determined by the number of negative terms. It is $n-k-1$ if $k<i$ and $n-k$ if $k>i+1$. The number of negative terms in $D_i^n$ is $n-i$.

Then the sign of $\frac{d}{dt}\tilde L_i^n(\tau_k)$ is determined by $2n-i-k-1$ or $2n-i-k$ negative terms with positive weights. Since $2n$ is always even, the sign is determined by $i$ and $k$.

Similarly, the sign of $\frac{d}{dt}\tilde L_{i+1}^n(\tau_k)$ is determined by $i$ and $k$ and is opposite to that of $\frac{d}{dt}\tilde L_i^n(\tau_k)$.

By the definition of $F_i^n$, the sign of $\frac{d^2}{dt^2}F_i^n(\tau_k)$, $k\ne i,k\ne i+1$, is as follows.

| Position of $k$ | $i+k$ | Sign of $\frac{d^2}{dt^2}F_i^n(\tau_k)$ |
| --- | --- | --- |
| $k<i$ | even | $(-1)^{i+k+1}=-1$ |
| $k<i$ | odd | $(-1)^{i+k+1}=+1$ |
| $k>i+1$ | even | $(-1)^{i+k}=+1$ |
| $k>i+1$ | odd | $(-1)^{i+k}=-1$ |

## Endpoint indices $i=0,n$

The endpoint definitions $G_0^n=\frac12$ and $G_{n+1}^n=-\frac12$ are constant. Hence

$$
(F_0^n)''(t)=-(\tilde L_1^n)'(t),\qquad
(F_n^n)''(t)=(\tilde L_n^n)'(t).
$$

For $k=2,\dots,n$, the numerator of $(\tilde L_1^n)'(\tau_k)$ has $n-k$ negative factors, while $D_1^n$ has $n-1$. Since $\omega_1^n>0$,

$$
\operatorname{sgn}(F_0^n)''(\tau_k)=(-1)^k.
$$

For $k=1,\dots,n-1$, the numerator of $(\tilde L_n^n)'(\tau_k)$ has $n-k-1$ negative factors, while $D_n^n$ has none. Since $\omega_n^n>0$,

$$
\operatorname{sgn}(F_n^n)''(\tau_k)=(-1)^{n-k-1}=(-1)^{n+k+1}.
$$

Thus the table above also holds at $i=0$ (only $k>i+1$) and $i=n$ (only $k<i$). Together, it covers $i=0,\dots,n$ and the applicable node critical points.

<br><br><br>

# Position of $t^*$

- t* is defined under $\omega \neq 0$
- $D_i^n \neq 0$ with the distinct node
- So, $\alpha \neq 0$ and $\beta \neq 0$
- $
\alpha\beta=\frac{\omega_{i+1}^n}{D_i^n}\frac{\omega_i^n}{D_{i+1}^n}.
$


$t^*$ divides the segment from $\tau_i$ to $\tau_{i+1}$ in the ratio $\alpha:-\beta$. We can consider 3 cases.
1. $t^*$ lies on $\tau_i$ or $\tau_{i+1}$.
2. $t^*$ is the internal division point of  $\tau_i$ and $\tau_{i+1}$.
3. $t^*$ is the external division point of  $\tau_i$ and $\tau_{i+1}$. 

### 1. $t^*$ lies on $\tau_i$ or $\tau_{i+1}$ 
Since $\alpha \neq 0$ and $\beta \neq 0$, this case is impossible.

### 2. $t^*$ is the internal division point of  $\tau_i$ and $\tau_{i+1}$.
$
   \iff \alpha(-\beta) > 0 
\\ \iff \alpha\beta < 0
\\ \iff \omega_i^n \omega_{i+1}^n > 0 \qquad \because D_i^n D_{i+1}^n < 0
\\ \iff \omega_i^n$ and $ \omega_{i+1}^n$ have the same sign.

### 3. $t^*$ is the external division point of  $\tau_i$ and $\tau_{i+1}$.
$
   \iff \alpha(-\beta) < 0 
\\ \iff \alpha\beta > 0
\\ \iff \omega_i^n \omega_{i+1}^n < 0 \qquad \because D_i^n D_{i+1}^n < 0
\\ \iff \omega_i^n$ and $ \omega_{i+1}^n$ have different signs.
$\\ \Longrightarrow \exist \omega <0$ 

### useful facts
1. Under the positive weights condition, $t^*$ is always the internal division point
2. With distinct nodes, if $t^*$ is the external division point, then there exists a nonpositive weight.

<br><br><br>

# New notation of $t^*$

Under the positive weight condition, $t^*$ is the internal point of $\tau_i$ and $\tau_{i+1}$. So from now on, we denote $t^*$ by $\tau_{i+0.5}$, because in many examples it is more understandable. It does not mean the geometrical midpoint of $\tau_i$ and $\tau_{i+1}$.

<br><br><br>

# $t^*$ is the local maximizer

For $i=1,\dots,n-1$, assume that all quadrature weights are positive. By the interpolation property of the Lagrange polynomials,

$$
\frac{d}{dt}F_i^n(\tau_i)=\frac{1}{\omega_i^n}>0,
\qquad
\frac{d}{dt}F_i^n(\tau_{i+1})=-\frac{1}{\omega_{i+1}^n}<0.
$$

Since $t^*=\tau_{i+0.5}$ is the unique critical point in $(\tau_i,\tau_{i+1})$, the derivative changes from positive to negative at $t^*$. Therefore, $t^*$ is a local maximum point of $F_i^n$.

<br><br><br>

# $t^*$ is not necessarily a global maximizer

Even if the nodes are symmetric about zero, include both endpoints, and have positive quadrature weights, $t^*$ need not be a global maximizer on $[-1,1]$.

For example, consider the six nodes

$$
(\tau_1,\dots,\tau_6)=\left(-1,-\frac12,-\frac15,\frac15,\frac12,1\right).
$$

Their quadrature weights are all positive:

$$
(\omega_1^6,\dots,\omega_6^6)=\left(\frac{17}{108},\frac{128}{189},\frac{125}{756},\frac{125}{756},\frac{128}{189},\frac{17}{108}\right).
$$

The corresponding basis polynomial is

$$
F_3^6(t)=-\frac{25}{8}(t^2-1)^2(8t^2+1).
$$

Its critical point between $\tau_3$ and $\tau_4$ is $t^*=0$, which is a local maximum. However,

$$
F_3^6(0)=-\frac{25}{8}<0=F_3^6(-1)=F_3^6(1).
$$

Thus $t^*$ is not a global maximizer on $[-1,1]$.

![Basis functions for the symmetric six-node counterexample](critical_point_counterexample.png)


<br><br><br>

<h1 style="background-color: rgb(134, 227, 182); color: rgb(0, 0, 0)"> Minimize Variance Problem </h1>

# Minimizing the variance of the maxima of F

We make the following assumptions.
- If the maxima of the basis functions are close to one another, the curve is more predictable.
- Likewise, the smaller the variance of the maxima, the more predictable the curve.
- If the maxima of the inner basis functions are close to one another, their minima are also close.


<br><br><br>

# Definitions

### Inner basis functions
We call $F_1^n, F_2^n, \dots, F_{n-1}^n$, i.e., the basis functions excluding the endpoint indices $i=0$ and $i=n$, the *inner basis functions*.

### Maximum of a basis function
Let $\tau = (\tau_1, \dots, \tau_n)$ be a strictly increasing sequence of nodes. <br>
For $i = 0, 1, \dots, n$, define the maximum of $F_i^n$ by
$$M_i^n(\tau) = \max_{t\in[-1,1]} F_i^n(t).$$

### Average and variance of the maxima
Define the average $\bar{M}(\tau)$ and the variance $\operatorname{Var}(\tau)$ of the maxima of the inner basis functions by

$$
\bar M(\tau)=\frac{1}{n-1}\sum_{i=1}^{n-1}M_i^n(\tau),\qquad
\operatorname{Var}(\tau)=\frac{1}{n-1}\sum_{i=1}^{n-1}\bigl(M_i^n(\tau)-\bar M(\tau)\bigr)^2 .
$$


<br><br><br>

# 탐색 방법

### 고려 사항
- 목적함수(분산)은 $\tau$에 따라 연속이지만 미분 불가능
- 즉, Gradient Discnet 사용 어려움 
- 분산을 하나 평가하는데 비용은 100ms 미만

### 후보1 : grid search
자유변수들의 격자점을 모두 탐색

장점: 국소최솟값 문제를 고민 안해도 됨
단점: 자유변수 개수 n에 따라 $격자 개수^n $만큼 탐색하여 차수가 증가할수록 비용이 많이 듦

### 후보2 : compass search
노드를 자유변수 하나만 +h, -h 로 움직여 후보를 생성
자유 변수 n에 따라 총 2n개의 후보를 살펴보고 분산이 최소가 되는 노드로 이동

장점: 간단한 구현, 낮은 비용
단점: 대각선에 낮은 비용이 위치하는 경우, 놓칠 수 있음

### 후보3 : compass search with diagonal
compass search를 대각선을 포함하여 실행
자유 변수 n개에 따라 총 $3^n - 1$개의 후보를 살펴봄

장점: 간단한 구현
단점: compass search보다 비싼 비용

### 후보4 : Mesh Adaptive Direct Search
n차원 격자점의 집합(mesh)의 간격과 poll의 크기, mesh의 후보점들을 기준에 따라 정하여 탐색하는 방법

- 공부 더 필요함

<br><br><br>

# Compass Search

1. 차수 n, 탐색 간격 h, 종료 조건 h_limit, 초기 노드 initial node 를 입력
2. 초기 노드의 분산을 계산
3. 하나의 자유변수에 +h, -h 연산을 수행한 후보 노드에 대해 분산을 계산
4. 3을 모든 자유변수에 대해 반복
5. 초기 분산과 후보 분산들 중 최솟값 탐색
6. 최솟값이 나오는 노드가 초기 노드일 경우, 탐색 간격을 10분의 1로 조정
7. 최솟값이 나오는 노드가 탐색 노드일 경우, 해당 노드를 초기 노드로 설정

![alt text](image.png)

---

### n = 4

![alt text](image-1.png)

  93    0.000000000100000 2.333690955214700e-19  [-1.00000000000000000000, -0.42818182699999973728, 0.42818182699999973728, 1.00000000000000000000]
[case3] a=0.42818183, sympy check: 2.333688426932632e-19

---

### n=5

![alt text](image-2.png)
[ESC] a=0.61525633, sympy check: 2.467236153920099e-19
[case1] a=0.36343098, sympy check: 1.971465173421264e-18
[case2] a=0.36343098, sympy check: 1.971466108732261e-18
[case3] a=0.61525634, sympy check: 4.025324281057727e-21
[case4] a=0.61525634, sympy check: 4.025370066276577e-21
[case5] a=0.82468554, sympy check: 8.877902265598494e-18


![alt text](image-3.png)
![alt text](image-4.png)
![alt text](image-5.png)

---

### n=6

![alt text](image-6.png)
![alt text](image-7.png)
![alt text](image-8.png)

sympy check: 1.806596360179689e+00
sympy check: 1.356292541170561e-02
sympy check: 3.425611140052588e-08
sympy check: 1.591464508424902e-10

![alt text](image-9.png)
![alt text](image-10.png)

---

### 관찰 결과

분산의 국소 최솟값이 여러개 존재하더라도 전역 최솟값은 양수 영역에서 나타날 것이다.
-> n = 7 인 경우, 그렇지 않음이 나타남
-> 1. compass serach의 한계로 인해 명제가 참이지만 그렇게 보이지 않음
-> 2. 전역 최솟값은 여러개 존재할 수 있음

그렇다면 새로운 추측
양수 가중치 영역에서 국소 최솟값은 하나이다?

- 가중치가 양수인 공간은 연결 공간이다.
- 가중치가 양수인 공간은 하나의 국소 최솟값을 갖는다.
- 가중치가 양수인 공간의 국소 분산 최솟값은 t*에 의해서 발생한다.
- t* 에서 항상 최댓값이 결정되는 공간은?
- 분산의 국소 최솟값은 0이다.
- compass search의 오차는 어느정도인가? 얼마나 신뢰해야하는가 → 아래 "Compass Search의 신뢰도" 참고
- 양수 가중치에서 갖는 최댓값이 차수와 무관하게 0.87 $\cdots$
GL basis의 상한과 비슷?
![alt text](image-11.png)
![alt text](image-12.png)
<br><br><br>

# Compass Search의 신뢰도 (2026-11-13, n = 6)

### 실험 방법
- 노드 $(-1,-b,-a,a,b,1)$, 기존 설정 그대로 `compass_search(6, h=0.001, h_limit=1e-6)`
- 시작점: $0.05$ 간격 격자에서 $a<b$인 171개
- 정답지: `variance_cache/n6/N1000`의 Fraction 정확값 격자 (간격 0.001)
- 각 도착점에서 반지름 $10^{-3}, 10^{-4}, 10^{-5}$ 원 위 720방향의 값을 비교해, 더 낮은 값이 있으면 "가짜 수렴"으로 판정
- 도착점에서 16방향 pattern search(간격을 절반씩 줄임)로 다시 내려가 실제 극소점과 비교

### 도착점 분류 (171회)

| 도착점 $(a, b)$ | 횟수 | compass가 보고한 분산 | 실제 극소값 | 가중치 |
| --- | --- | --- | --- | --- |
| $(0.25408, 0.72590)$ | 63 (37%) | $1.4\times10^{-10}\sim1.8\times10^{-10}$ | $7\times10^{-27}$ | 양수 |
| $(0.32604, 0.87521)$ | 25 (15%) | $2.7\times10^{-8}\sim3.4\times10^{-8}$ | $7\times10^{-25}$ | 음수 존재 |
| $(0.46247, 0.99560)$ | 2 (1%) | $3.5\times10^{-9}$ | $2\times10^{-24}$ | 음수 존재 |
| $a \approx b$ (노드 충돌) | 81 (47%) | $1.4\times10^{-2}\sim6.4$ | 정의 안 됨 | - |

### 관찰 결과
1. **수치 오차는 문제가 아니다.** numpy 분산과 Fraction 정확값의 차이는 최대 $4\times10^{-14}$ (상대오차 최대 $10^{-10}$).
2. **노드 충돌 버그 (47%).** $b-a \approx 10^{-16}$인 상태로 끝난 경우가 81회.
   `n_ >= node[right + 1]` 검사에서 $a+h$가 반올림 때문에 $b$보다 $10^{-16}$만큼 작게 계산되면 통과해 버린다.
   이 점들은 서로 다른 노드라는 가정을 깨므로 결과로 쓰면 안 된다. (자세한 내용은 아래 "노드 충돌 버그" 참고)
3. **극소점 근처에서 멈춘다 (가짜 수렴).** 세 내부 극소점 모두 반지름 $10^{-5}$ 안에 compass가 보고한 값보다 낮은 점이 있었다.
   대각선 방향 탐색으로 계속 내려가면 분산이 $10^{-25}$ 수준(부동소수점 한계)까지 떨어진다.
   → n=6에서 찾은 세 국소 최솟값은 모두 **0**으로 보인다 ("분산의 국소 최솟값은 0이다" 추측을 지지).
   → compass가 보고하는 $10^{-10}, 10^{-8}$ 같은 값은 실제 최솟값이 아니고, 좌표축 방향으로 더 못 내려가는 지점의 값일 뿐이다.
4. **위치 정확도는 $10^{-4}\sim10^{-5}$ 수준이다.** `h_limit = 1e-6`이어도 도착점과 실제 극소점의 거리는 $3\times10^{-6}\sim9\times10^{-5}$.
   예) case3 optimize $(0.326110, 0.875282)$ → 실제 $(0.326037, 0.875210)$. 소수점 5자리 이후는 의미가 없다.
5. **어느 극소점에 도착할지는 시작점에 달려 있다.** 양수 가중치 극소점에 도착한 것은 37%뿐이다.
6. **비용.** 반복 횟수 중앙값 544, 최대 21122 (1회 실행 중앙값 약 11초).

### 결론: 얼마나 믿을 수 있는가
- 어느 **분지(basin)** 에 있는지, 극소점의 대략적 위치($10^{-4}$)는 믿을 수 있다.
- 최솟값 자체와 소수점 5자리 이후의 위치는 믿을 수 없다.
- $a\approx b$로 끝난 결과는 버려야 한다.

### 개선 방향
- 경계 검사에 여유를 둔다 (예: `n_ >= node[right + 1] - eps`).
- 좌표축 방향만이 아니라 대각선 / 여러 방향을 함께 본다 (후보3, 후보4 MADS).
- 최종 점만 Nelder–Mead나 여러 방향 pattern search로 다듬는다.

<br><br><br>

# 노드 충돌 버그

### 재현 (시작점 $(a, b) = (0.10, 0.15)$)

| step | $a$ | $b$ | $b-a$ | 분산 | 다음 $h$ |
| --- | --- | --- | --- | --- | --- |
| 49 | 0.14900000000000005 | 0.15 | $1.0\times10^{-3}$ | 3.0012 | 0.001 |
| 51 | 0.14910000000000004 | 0.15 | $9.0\times10^{-4}$ | 2.9971 | 0.0001 |
| $\vdots$ | | | | | |
| 59 | 0.14989999999999995 | 0.15 | $1.0\times10^{-4}$ | 2.9647 | 0.0001 |
| 60 | 0.14999999999999994 | 0.15 | $5.6\times10^{-17}$ | 2.9606 | 0.0001 |

이후 더 움직이지 못하고 $h$만 줄어들다 종료한다.

### 원인: 경계 검사의 반올림
```python
n_ = node[right] + sign*h
if n_ <= 0 or n_ <= node[right - 1]: continue
if n_ >= node[right + 1]: continue   # a + h >= b 이면 건너뜀
```
- 0.001, 0.0001은 2진수로 정확히 표현되지 않아 `a += h`를 반복하면 마지막 자리 오차가 쌓인다.
- step 59에서 수학적으로는 $a + 10^{-4} = b$이므로 검사에 걸려야 한다.
- 하지만 float로는 $a + 10^{-4} = 0.14999999999999994 < b$가 되어 통과한다.
- 결과적으로 두 노드가 약 $10^{-16}$ 간격으로 거의 겹친다.

### 왜 그 자리에 갇히는가
- 이 시작점에서는 $a=b$ 쪽으로 갈수록 분산이 줄어든다. 즉 정의역의 벽 $a=b$ 쪽이 내리막이다.
- 붙고 나면 $a$를 올리는 이동과 $b$를 내리는 이동은 검사에 걸리고, $a$를 내리거나 $b$를 올리면 분산이 커진다.
- 그래서 $h$만 줄어들다 종료하고, 이 점을 극소점으로 보고한다.

### 노드가 가까워지면 생기는 일
- 간격을 $d$라 하면 $D_i \propto d$이고, 가중치와 $L_i$는 $1/d$에 비례해 커진다.
  예) 위 도착점의 가중치: $(0.19,\ 7.3\times10^{15},\ -7.3\times10^{15},\ -7.3\times10^{15},\ 7.3\times10^{15},\ 0.19)$
- 큰 값들이 상쇄되어 $F_i$가 만들어지므로 float 계산의 상대오차는 대략 $10^{-16}/d$이다.
  - $d = 10^{-4}$이면 약 $10^{-12}$로 문제없다.
  - $d \approx 10^{-16}$이면 결과를 믿을 수 없다.
- 수학적으로는 분산이 유한한 극한으로 수렴한다 (간격 $10^{-3}, 10^{-6}, 10^{-9}, 10^{-12}$에서 2.9648 → 2.9606 → 2.9606 → 2.9606).
  즉 float 계산의 문제이고, Fraction으로 계산하면 정확하다.

### 같은 문제가 생길 수 있는 곳 (n = 6)

| 가까워지는 노드 | 조건 | 관련 검사 |
| --- | --- | --- |
| $a$와 $b$ | $b-a \to 0$ | `n_ >= node[right+1]`, `n_ <= node[right-1]` |
| $-a$와 $a$ | $a \to 0$ | `n_ <= 0` |
| $b$와 $1$ | $b \to 1$ | `n_ >= node[right+1]` |

171회 실험에서는 $a$–$b$ 충돌만 나타났다. n ≥ 7이면 $a$와 가운데 노드 0처럼 인접한 노드 쌍이 더 생긴다.

### 주의점
- **최솟값 결론에는 영향이 없다.** 충돌로 끝난 81회의 분산은 $1.4\times10^{-2}\sim6.4$로, 내부 극소점(약 0)과 자릿수부터 다르다.
  위치도 대각선 $a \approx b$ 위라서 양수 가중치 극소점 $(0.254, 0.726)$과 멀다. 분산값만 봐도 걸러진다.
- **시작점을 고를 때:** 실행의 47%가 이 벽으로 빠진다. $a$와 $b$가 충분히 떨어진 시작점(예: $b-a > 0.1$)을 쓰면 계산 시간을 아낄 수 있다.
  "양수 가중치 극소점 도달 37%"는 이 47%를 포함한 비율이다.
- **정답 격자가 없는 n에서:** n = 8, 9 이상에서는 경계 쪽 극소값이 얼마나 낮은지 아직 모른다.
  결과를 볼 때마다 `min(np.diff(node)) < 1e-6`인지 확인한다.
- **버그를 고쳐도 남는 문제:** eps 여유를 두면 충돌은 막지만, 이 시작점들은 벽 바로 앞(간격 약 $h$)에서 멈춘다.
  "경계에서 멈춤"으로 표시하고 버리거나 다른 시작점으로 다시 돌려야 한다.

### 수정 방법
```python
eps = 1e-12
if n_ <= eps or n_ <= node[right - 1] + eps: continue
if n_ >= node[right + 1] - eps: continue
```
종료 후 `min(np.diff(node)) < 2*h_limit`이면 경계 수렴으로 표시한다.

<br><br><br>



# 고민

1. 최댓값을 갖는 점이 여러개인 경우에는 어떤 점을 해당 $$F_i^n$$의 maximizer로 선택할 것인가 -> 근데 이걸 고민해야 되나?
2. 가중치가 양수임을 가정하면 내분점임을 보장할 수 있지만, 내분점과 외분점으로는 가중치를 보장할 수 없음. -> 이걸 잘 알아보면 양수 가중치 조건을 찾는데 도움이 될지도