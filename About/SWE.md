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

Найдем собственные числа.

$$
\det \begin{pmatrix}
u -\lambda & h \\
g & u - \lambda
\end{pmatrix} \text{=} (u - \lambda)^2 - g h \text{=} 0
$$

$$
\lambda_{1,2} = u \pm c
$$

где 

$$
c = \sqrt{g h}
$$

Тогда легко найти собственные векторы:

$$
\vec l_1 = 
\begin{pmatrix}
-c \\
h
\end{pmatrix}
$$

$$
\vec l_2 = 
\begin{pmatrix}
c \\
h
\end{pmatrix}
$$

Теперь мы умножаем исходную систему на собственные векторы слева.
