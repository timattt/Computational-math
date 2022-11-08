# Численное решение задач по квантмеху

### Условие

Найти ВКБ спектр частицы в потенциале

$$
U(x) = U_0 (\frac{a}{x} - \frac{x}{a})^2
$$

### Решение

![image](https://user-images.githubusercontent.com/25401699/200623173-28b0d3da-e717-4e6b-a43f-0a493f34a206.png)

![image](https://user-images.githubusercontent.com/25401699/200623208-d1288d79-8c0c-4c9f-9b07-298715170d17.png)

![image](https://user-images.githubusercontent.com/25401699/200623257-2163283d-e2f0-444e-85bc-8ec4aec97f96.png)

![image](https://user-images.githubusercontent.com/25401699/200623301-e1dff6b3-1f8c-4b51-ac11-1fb5e5417731.png)

### Численное построение волновых функций

Теперь, когда мы умеем находить энергии, мы легко можем численно проинтегрировать в трех областях.
В разрешенной:

$$
\psi(x) = \frac{C}{\sqrt{p(x)}}\sin (\frac{1}{\hbar}\int_{a}^{x}p(x')dx' + \frac{\pi}{4})
$$

В запрешенной слева:

$$
\psi(x) = \frac{C}{2\sqrt{|p(x)|}}\exp (-\frac{1}{\hbar}\int_{x}^{a} |p(x')| dx')
$$

В запрешенной справа:

$$
\psi(x) = \frac{C}{2\sqrt{|p(x)|}}\exp (-\frac{1}{\hbar}\int_{b}^{x} |p(x')| dx')
$$
