# Двумерное уравнение мелкой воды

Имеем вот такую систему.

$$
\begin{cases}
  \frac{\partial H}{\partial t} + \frac{\partial (H u)}{\partial x} + \frac{\partial (H v)}{\partial y} = 0 \\
  \frac{\partial (H u)}{\partial t} + \frac{\partial}{\partial x}(H u^2 + \frac{1}{2}gH^2) + \frac{\partial}{\partial y} (H u v) = 0 \\
  \frac{\partial (H v)}{\partial t} + \frac{\partial}{\partial x} (H u v) + \frac{\partial}{\partial y}(H v^2 + \frac{1}{2}gH^2) = 0
\end{cases}\
$$

Сделаем замену.

$$
\begin{pmatrix}
H \\
H u \\
H v
\end{pmatrix} =
\begin{pmatrix}
H \\
q \\
p
\end{pmatrix} = \vec U
$$

Тогда можно переписать в матричном виде.

$$
\frac{\partial \vec U}{\partial t} + 
\begin{pmatrix}
0 & 1 & 0 \\
g H - u^2 & 2u & 0 \\
-u v & v & u
\end{pmatrix}
\frac{\partial}{\partial x} \vec U
+
\begin{pmatrix}
0 & 0 & 1 \\
-u v & v & u \\
g H - v^2 & 0 & 2 v
\end{pmatrix}
\frac{\partial}{\partial y} \vec U = 0
$$
