# Уравнения мелкой воды

$$
 \frac{\partial h}{\partial t} + \frac{\partial (h u)}{\partial x} = 0
$$

$$
\frac{\partial (h u)}{\partial t} + \frac{\partial (h u^2)}{\partial x} + g h \frac{\partial h}{\partial x} = 0
$$


Перепишем в матричном виде:

$$
\begin{pmatrix}
\frac{\partial h}{\partial t} \\
\frac{\partial u}{\partial t}
\end{pmatrix}
+
\begin{pmatrix}
u & h \\
g & u
\end{pmatrix}
*
\frac{\partial}{\partial x}
\begin{pmatrix}
h \\
u
\end{pmatrix}
\text{=}
\vec 0
$$

<img src="https://user-images.githubusercontent.com/25401699/196695879-f8bc6d78-cf92-4ba6-80f1-7e2a0cf2b127.png" alt="drawing" width="500"/>


Найдем собственные числа.

![image](https://user-images.githubusercontent.com/25401699/196630089-df0ae447-6324-480b-ab75-341f91ce8a73.png)

Введем обозначение

$$
c = \sqrt{g h}
$$

Тогда легко найти собственные векторы:

![image](https://user-images.githubusercontent.com/25401699/196630369-6b198e9d-46f2-4864-ba0a-16c8d04cd565.png)

Теперь можем умножить на левые собственные векторы и ввести инварианты Римана.

$$
r = u + 2c
$$

$$
s = u - 2c
$$

Тогда имеем:

$$
\frac{\partial s}{\partial t} + (\frac{1}{4}r + \frac{3}{4}s)\frac{\partial s}{\partial x} = 0
$$

$$
\frac{\partial r}{\partial t} + (\frac{3}{4}r + \frac{1}{4}s)\frac{\partial r}{\partial x} = 0
$$
