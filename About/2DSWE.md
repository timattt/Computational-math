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

У нас появилось две матрицы.
Давайте решим спектральную задачу для обеих.

### Матрица при dx

$$
A = \begin{pmatrix}
0 & 1 & 0 \\
g H - u^2 & 2u & 0 \\
-u v & v & u
\end{pmatrix}
$$

Тогда собственные значения и их собственные векторы будут:

$$
\lambda_1 = u
$$

$$
\vec l_1 = \begin{pmatrix}
0 \\
0 \\
1
\end{pmatrix}
$$


$$
\lambda_2 = u - c
$$

$$
\vec l_2 = \begin{pmatrix}
1 \\
u - c \\
v
\end{pmatrix}
$$


$$
\lambda_3 = u + c
$$

$$
\vec l_3 = \begin{pmatrix}
1 \\
u + c \\
v
\end{pmatrix}
$$

### Матрица при dy

$$
B = \begin{pmatrix}
0 & 0 & 1 \\
-u v & v & u \\
g H - v^2 & 0 & 2 v
\end{pmatrix}
$$

Тогда собственные значения и их собственные векторы будут:

$$
\lambda_1 = v
$$

$$
\vec l_1 = \begin{pmatrix}
0 \\
1 \\
0
\end{pmatrix}
$$


$$
\lambda_2 = v - c
$$

$$
\vec l_2 = \begin{pmatrix}
1 \\
u \\
v - c
\end{pmatrix}
$$


$$
\lambda_3 = v + c
$$

$$
\vec l_3 = \begin{pmatrix}
1 \\
u \\
v + c
\end{pmatrix}
$$

### Численное решение

Будем решать нашу задачу в два этапа.
На первом этапе решаем задачу вот такую:

$$
\frac{\partial}{\partial t} \vec U + A \frac{\partial}{\partial x} \vec U = 0
$$

Будем умножать эта равенство последовательно на i-ый левый собственный вектор матрицы A.

$$
(\vec l_i * \frac{\partial}{\partial t} \vec U) + (\vec l_i * A \frac{\partial}{\partial x} \vec U) = 0
$$

Преобразуем исходя из определения собственного вектора.


$$
(\vec l_i * \frac{\partial}{\partial t} \vec U) + (\vec l_i * \lambda_i \frac{\partial}{\partial x} \vec U) = 0
$$

Теперь будем апроксимировать производную по пространству оператором из одномерного случая.

$$
(\vec l_i * \frac{\partial}{\partial t} \vec U) + (\vec l_i * \lambda_i L_{\lambda_i}[\vec U]) = 0
$$

Теперь апросимируем по времени.

$$
(\vec l_i * \frac{\hat{\vec U} - \vec U}{\tau) + (\vec l_i * \lambda_i L_{\lambda_i}[\vec U]) = 0
$$
