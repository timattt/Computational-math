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
(\vec l_i * \frac{\hat{\vec U} - \vec U}{\tau} + (\vec l_i * \lambda_i L_{\lambda_i}[\vec U]) = 0
$$

Теперь пусть мы уже посчитали сетку на k-ом временном слое. Пусть эта сетка U. Тогда отсюда мы легко найдем U с шляпкой.
Для этого запишем последнее уравнение для всех трех собственных векторов.

$$
\begin{cases}
   p_t + (\vec l_1 * \lambda_1 L_{\lambda_1}[\vec U]) = 0 \\
   H_t + q_t (u-c) + p_t v + (\vec l_2 * L_{\lambda_2}[\vec U] = 0 \\
   H_t + q_t (u+c) + p_t v + (\vec l_3 * L_{\lambda_3}[\vec U] = 0
\end{cases}\
$$

Эту систему можно решить и получить производные по времени в данный момент времени.
После чего, зная апроксимацию 

$$
\vec U_t = \frac{\hat{\vec U} - \vec U}{\tau}
$$

Мы получим значения p, q, H со шляпкой. То есть для нового временного слоя.
Мы нашли U со шляпкой для этой задачи. Теперь будем решать такую задачу:

$$
\frac{\partial}{\partial t} \vec U + B \frac{\partial}{\partial x} \vec U = 0
$$

Действуя аналогично предыдущему случаю, получим следующую систему для нахождения k+1 временного слоя.

$$
\begin{cases}
   q_t + (\vec l_1 * \lambda_1 L_{\lambda_1}[\vec U]) = 0 \\
   H_t + q_t u + p_t (v-c) + (\vec l_2 * L_{\lambda_2}[\vec U] = 0 \\
   H_t + q_t u + p_t (v+c) + (\vec l_3 * L_{\lambda_3}[\vec U] = 0
\end{cases}\
$$

Отсюда снова можем найти q, p, H для следующего временного слоя.

Теперь сделаем следующее:
1. У нас есть k-ый временной слой Uij.
2. Мы решим первую задачу на всех узлах сетки, используя этот слой и получим новый слой Vij.
3. Мы решим вторую задачу на всех узлах сетки, используя промежуточный слой Vij и получим итоговый новый слой в основной задаче.

То есть по сути своей мы как бы вспахиваем двумерную сетку xy как поле. Сначала, решая первой задачей, линиями параллельными оси x.
А потом, решая второй задачей, линией параллельными оси y.

С граничными условиями на каждой такой линии разбираемся аналогично одномерному случаю. То есть смотрим на лямбды и на то, будет ли схема явной или неявной.
Отметим, что в итоговых системах два последних уравнения помогают найти на границе H. А первое согласовать u и v с их граничными условиями, но
если мы на границе требуем равенство скорости нулю, то первое уравнение нам там не понадобится.

# Тестирование

Имеем сетку XoY 30 на 30.
И 200 вперед по времени.

$$
\begin{cases}
h = 0.01\\
\tau = 0.001\\
N = 30\\
M = 200\\
L = N * h = 0.3
\end{cases}
$$

## Экспонента в центре

### Начальные условия

$$
\begin{cases}
H(x, y) = 1 + 0.3*h* \exp(-500*((x-\frac{L}{2})^2 + (y-\frac{L}{2})^2)) \\
u(x, y) = v(x, y) = 0
\end{cases}
$$

![image](https://user-images.githubusercontent.com/25401699/199317900-36ad95fc-f74b-4832-a31f-3734516916e0.png)

### Граничные условия

На границе сетки:
$$
u = v = 0
$$

### Процесс

![test](https://user-images.githubusercontent.com/25401699/199316580-75828ff5-873b-49bb-8d44-f427816d98fa.gif)

