## Определение корней уравнения

### Метод двоичного поиска

Мы тут люди скромные и всякие эти ваши понтовые названия алгоритмов вроде "дихотомия" не используем!

Алгоритм такой: у нас есть функция **f** и отрезок **[a, b]** на этом отрезке у этой функции проживает по постоянной прописке корень. Делим отрезок много раз рекурсивно пополам и находим корень.

Подробности [тут](http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_%D0%B4%D0%B8%D1%85%D0%BE%D1%82%D0%BE%D0%BC%D0%B8%D0%B8)

### Метод Ньютона (метод секущих)

![](https://github.com/timattt/Project-computational-math/blob/master/Images/Newton_tangents.png)

**Определение красности корня:**

Пусть **p** - кратность корня, тогда:
```
p = 1 / (1 - Q)
Q = lim(Qn) where n -> oo
Qn = (Xn+1 - Xn) / (Xn - Xn-1)
```

### Метод Ньютона (метод секущей)

![](https://github.com/timattt/Project-computational-math/blob/master/Images/Newton_tangent.png)

### Метод Ньютона (метод хорд)

![](https://github.com/timattt/Project-computational-math/blob/master/Images/Newton_chords.png)

### Метод Ньютона (многомерный случай)

Имеем вектор-функцию **f**. Имеет ее матрицу якоби **F'**.
Тогда для поиска корней справедлива формула рекурсивная:

![](https://github.com/timattt/Project-computational-math/blob/master/Images/Newton_multi.png)

Где **Xk** - векторная последовательность.

Но, как можно заметить - это формула содержит обратную матрицу. Давайте получим другую формулу для вычисления **Xk+1**. Сведем задачу к решению СЛАУ.

```
F'(Xk)*Xk+1 = F'(Xk)*Xk - f(Xk)
A = F'(Xk)
B = F'(Xk)*Xk - f(Xk)

A*Xk = B
```

Теперь можно решить эту СЛАУ и получить нужное значение быстрее.
