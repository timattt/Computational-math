# Уравнения мелкой воды

$$
 \frac{\partial H}{\partial t} + \frac{\partial (H u)}{\partial x} = 0
$$

$$
\frac{\partial (H u)}{\partial t} + \frac{\partial (H u^2)}{\partial x} + g H \frac{\partial H}{\partial x} = 0
$$


Перепишем в матричном виде:

$$
\begin{pmatrix}
\frac{\partial H}{\partial t} \\
\frac{\partial u}{\partial t}
\end{pmatrix}
+
\begin{pmatrix}
u & H \\
g & u
\end{pmatrix}
*
\frac{\partial}{\partial x}
\begin{pmatrix}
H \\
u
\end{pmatrix}
\text{=}
\vec 0
$$

Найдем собственные числа.

$$
\det \begin{pmatrix}
u -\lambda & H \\
g & u - \lambda
\end{pmatrix} \text{=} (u - \lambda)^2 - g H \text{=} 0
$$

$$
\lambda_{1,2} = u \pm c
$$

где 

$$
c = \sqrt{g H}
$$

Тогда легко найти собственные векторы:

$$
\vec l_1 = 
\begin{pmatrix}
-c & H
\end{pmatrix}
$$

$$
\vec l_2 = 
\begin{pmatrix}
c & H
\end{pmatrix}
$$

Теперь мы умножаем исходную систему в матричном виде на собственные векторы слева. Произведение матрицы и ее собственного вектора преобразуем
по определению собственного вектора. И получаем два выражения.

$$
-(H_t + \lambda_1 H_x) c + (u_t + \lambda_1 u_x) H = 0
$$

$$
(H_t + \lambda_2 H_x) c + (u_t + \lambda_2 u_x) H = 0
$$

Теперь во имя устойчивости при апроксимации по координате введем такой оператор:

$$
L_{\lambda_k}[U]_i^{j} = [\lambda_i^j]_{+}\frac{U_i^j - U_{i-1}^j}{h} + [\lambda_{i}^{j}]_{-} \frac{U_{i-1}^j - U_{i}^{j}}{h}
$$

И будем аппроксимировать таким образом:

$$
\lambda_k U_x \approx L_{\lambda_k}[U]_i^{j}
$$

$$
U_t \approx \frac{U_i^{j+1} - U_i^j}{\tau}
$$


