# Уравнения мелкой воды

## Численное решение

Имеем вот такую систему.

$$
\begin{cases}
   \frac{\partial H}{\partial t} + \frac{\partial (H u)}{\partial x} = 0\\
  \frac{\partial (H u)}{\partial t} + \frac{\partial (H u^2)}{\partial x} + g H \frac{\partial H}{\partial x} = 0
\end{cases}\
$$

Где u - скорость слоя жидкости. H - высота слоя жидкости.

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
L_{\lambda_k}[U]_i^{j} = [\lambda_i^j]_+ \frac{U_i^j - U_{i-1}^j}{h} + [\lambda_{i}^{j}]_- \frac{U_{i-1}^j - U_{i}^{j}}{h}
$$

где

$$
[\lambda_i^j]_+ = \frac{\lambda_i^j + |\lambda_i^j|}{2}
$$

$$
[\lambda_i^j]_- = \frac{\lambda_i^j - |\lambda_i^j|}{2}
$$

Отметим, определение устойчивости слудует из спектрального критерия Неймана.
Будем аппроксимировать таким образом:

$$
\lambda_k U_x \approx L_{\lambda_k}[U]_i^{j}
$$

$$
U_t \approx \frac{U_i^{j+1} - U_i^j}{\tau}
$$

Подставляем наши авпроксимации и получаем систему из двух неизвестных.

$$
H_i^{j+1} * [-\frac{c_i^j}{\tau}] + u_i^{j+1} * [\frac{H_i^j}{\tau}] = -\frac{c_i^j H_i^j}{\tau} + c_i^{j} L_{\lambda_1}[H]_i^{j} + \frac{H_i^j u_i^j}{\tau} - H_i^j L_{\lambda_1}[u]_i^{j}
$$

$$
H_i^{j+1} * [\frac{c_i^j}{\tau}] + u_i^{j+1} * [\frac{H_i^j}{\tau}] = \frac{c_i^j H_i^j}{\tau} - c_i^{j} L_{\lambda_2}[H]_i^{j} + \frac{h_i^j u_i^j}{\tau} - h_i^j L_{\lambda_2}[u]_i^{j}
$$

Имеем вот такую сетку.

![image](https://user-images.githubusercontent.com/25401699/197346542-cf91ded9-06f2-43d7-9b17-bbdce2298b6e.png)


Здесь зеленым обозначены начальные условия.
Внутренние точки, которые на схеме обозначены просто точками мы можем найти используя предыдущую систему.
Однако, что делать на крайних точках. Из формулы следует, что мы можем выйти за пределы сетки.
Введем граничные условия.

$$
u(x = 0) = u(x = L) = 0
$$

Тогда нам будет достаточно только одного уравнения из системы.
Из анализа системы очевидно следует, что при $$\lambda_k > 0$$ будет явный уголок. То есть берется точка $$i-1$$
Соответственно на левом конце выбираем уравнение, где $$\lambda_k < 0$$ А на правом то, где $$\lambda_k > 0$$

Если на краю нужного уравнения не нашлось, значит получается, что $$|u| > c$$
то есть скорость движения жидкости больше скорости звука, чего быть не может. Поэтому нужно выбирать данные для задачи с умом. С учетом условий УМВ.

## Тестирование

### Колокольная капля

$$L = 1$$

$$u(t = 0) = 0$$

$$H(x, t = 0) = 1 + e^{-500(x-\frac{L}{2})^2}$$

$$u(x = 0, t) = u(x = L, t) = 0$$

![Recording 2022-10-23 at 12 13 26](https://user-images.githubusercontent.com/25401699/197384116-8a98acf4-2a67-475b-af4e-72123cc7eec7.gif)

### Равномерный синус

$$L = 1$$

$$u(t = 0) = 0$$

$$H(x, t = 0) = 0.1 + 0.01\sin(10*\frac{x}{L})$$

$$u(x = 0, t) = u(x = L, t) = 0$$

![Recording 2022-10-23 at 12 25 43](https://user-images.githubusercontent.com/25401699/197384494-035cf517-f02f-498e-9e8c-4af38c5496a7.gif)

Еще вот есть вариант без анимации:

![image](https://user-images.githubusercontent.com/25401699/197385064-111cfd39-23c0-4c81-8190-7a7245933092.png)

## Аналитическое решение

### Линеаризация

Пусть 

$$
u(x, t) = u_0 + \delta u
$$

$$
H(x, t) = H_0 + \delta H
$$

положим тут

$$
H_0 = 1
$$

$$
u_0 = 0
$$

Теперь подставим это в нашу исходную систему.

$$
\begin{cases}
   \frac{\partial}{\partial t}(\delta H) + \frac{\partial}{\partial x}(\delta u (1 + \delta H)) = 0\\
  \frac{\partial}{\partial t}(\delta u + \delta u \delta H) + \frac{\partial}{\partial x}((1 + \delta H)\delta u^2) + g (1+\delta H) \frac{\partial}{\partial x}(\delta H) = 0
\end{cases}\
$$

Или в матричном виде:

$$
\begin{pmatrix}
\frac{\partial \delta H}{\partial t} \\
\frac{\partial \delta u}{\partial t}
\end{pmatrix}
+
\begin{pmatrix}
0 & 1 \\
g & 0
\end{pmatrix}
*
\frac{\partial}{\partial x}
\begin{pmatrix}
\delta H \\
\delta u
\end{pmatrix}
\text{=}
\vec 0
$$

Решим спектральную задачу.

$$
\lambda_1 = \sqrt{g}
$$

$$
\vec l_1 = 
\begin{pmatrix}
1 \\
\sqrt{g}
\end{pmatrix}
$$

$$
\lambda_2 = -\sqrt{g}
$$

$$
\vec l_2 = 
\begin{pmatrix}
-1 \\
\sqrt{g}
\end{pmatrix}
$$
