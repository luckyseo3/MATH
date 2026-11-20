# 수치 미분

## 전방차분
$$
f(x+h) = f(x) + hf'(x) + \frac{h^2}{2}f''(x) + \frac{h^3}{6}f'''(y) \\

\frac{f(x+h)-f(x)}{h} = f'(x) + \frac{h}{2}f''(x) + \frac{h^2}{6}f'''(y) = f'(x) + O(h) + o(h) = f'(x) + O(h)
$$

빅오 함수 O(h)와 스몰오 함수 o(h)의 차이는 다음과 같다.

O(h) : $\lim_{h \to 0}\frac{O(h)}{h} = k$ for $ k \in \mathbb{R}$  
o(h) : $\lim_{h \to 0}\frac{O(h)}{h} = 0$

k = 0 일 수도 있기 때문에 o(h) 함수는 O(h) 함수이기도 하다. 

## 후방 차분
$$
f(x-h) = f(x) - hf'(x) + \frac{h^2}{2}f''(x) - \frac{h^3}{6}f'''(y) \\

\frac{f(x)-f(x-h)}{h} = f'(x) - \frac{h}{2}f''(x) + \frac{h^2}{6}f'''(y) = f'(x) + O(h) + o(h) = f'(x) + O(h)
$$


## 중앙 차분
$$
f(x+h)-f(x-h) = 2hf'(x) + 2\frac{h^3}{6}f'''(y)\\
\frac{f(x+h)-f(x-h)}{2h} = f'(x)  + \frac{h^2}{6}f'''(y) = f'(x) + O(h^2)
$$

정리하면 전방 차분과 후방 차분은 $O(h)$ 수준이고,
중앙 차분은 $O(h^2)$ 수준이다.















<br><br><br><br><br>
# 절단오차

절단오차는 함수를 n차 테일러 급수로 근사하며 발생하는 오차를 말한다.

## 테일러 정리
(a,b) 에서 $C^{n+1}$인 함수 $f$에 대해서 
(a,b)의 중심 c에서 전개한 n차 테일러 다항식과 함수 $f$의 차이는 $\frac{f^{(n+1)(y)}}{(n+1)!}(x-c)^{(n+1)}$ for $y\in(c, x)$ 이다.

이때 $\frac{f^{(n+1)(y)}}{(n+1)!}(x-c)^{(n+1)} \leq \operatorname{inf}_{t\in(a,b)}\frac{f^{(n+1)(t)}}{(n+1)!}(x-c)^{(n+1)} \lt \infin$ 일 경우,
원함수와 근사식의 오차에는 상한이 존재한다.

