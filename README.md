# Избранные задачи вычислительной математики
В этом репозитории представленны избранные задачи вычислитеьной математики.

## Численное дифференцирование
Задача такая - есть непрерывная функция, нужно взять ее первую и вторую производные. Двумя способами, так чтобы точность была пропрциональна разбиению и его квадрату.

### Простейший случай
Самый простой вариант - это первая производная с первым, где погрешность пропорциональна порядку разбиения.
```
f'n = [f(Xn+1) - f(Xn)] / h;
Xn - сетка аргумента.
f - наша непрерывная функция.
h - размер разбиения.
f'n - сетка производной функции.
```
Такая формула не будет работать на правом краю, поэтому там просто скопируем значение f'n-1.

### Большая точность
Если мы хотим большей точности, то из ряда тейлора можно вывести формулу, которая даст большую точность.
Опустим разложения в ряды. Вот [тут](https://ru.m.wikipedia.org/wiki/%D0%9A%D0%BE%D1%8D%D1%84%D1%84%D0%B8%D1%86%D0%B8%D0%B5%D0%BD%D1%82%D1%8B_%D1%84%D0%BE%D1%80%D0%BC%D1%83%D0%BB_%D1%87%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D0%BE%D0%B3%D0%BE_%D0%B4%D0%B8%D1%84%D1%84%D0%B5%D1%80%D0%B5%D0%BD%D1%86%D0%B8%D1%80%D0%BE%D0%B2%D0%B0%D0%BD%D0%B8%D1%8F) есть коэффициенты для более точных формул.

В этой статье есть три таблицы - коэффициенты вперед, назад и симметрично.
Если мы хотим найти f'n, то вперед будет обращаться к Xk только таким, у которых k >= n. Назад, соответственно у которых k <= n. А симметричный вариант использует равное колличество коэффициентов спереди и сзади.

Теперь давайте запишем формулы для первой производной со вторым порядком точности.   
***Коэффициенты вперед:***
```
f'n = [-3/2 * f(Xn) + 2 * f(Xn+1) - 1/2*f(Xn+2)] / h;
```
***Коэффициенты назад:***
```
f'n = [3/2 * f(Xn) + -2 * f(Xn-1) + 1/2*f(Xn-2)] / h;
```
Теперь можно спокойно считать производную.
Будем использовать формулу 'вперед' для всех, кроме двух правых значений. Для них используем формулу 'назад'.
Вот небольшой пример для
```
f(x) = sin(x)
```
На графике можно видеть несколько порядков точности производной. Видно, что второй порядок лучше, чем первый.

![](https://github.com/timattt/Project-computational-math/blob/master/Images/DiffExample.png)

### Эксперимент с порядком точности
Теперь давайте напишем программу, которая для каждой из вышеуказнных формул построит график зависимости погрешности от мелкости разбиения.
Погрешность сетки будем определять как максимум из всех отклонений значений сетки от значения функции в точках сетки.

Добавим еще в наш эксперимент два порядка вторых производных.   
Имеем четыре функции вида:
```
sigma(h) = h ^ p;
sigma - погрешность.
h - мелкость разбиения.
p - порядок.
```
Давайте экспериментально проверим, что для первых двух порядков точности первых двух производных. Точности соответсвенно p = 1 для первого порядка и p = 2 для второго порядка.   

Построим график:

![](https://github.com/timattt/Project-computational-math/blob/master/Images/GraphDiffRaw.png)

По оси X - мелкость разбиения. По оси Y - погрешность.
Можно заметить, что у разной точности разный изгиб линии.
Чтобы найти p, прологарифмируем обе стороны:
```
sigma(h) = h ^ p;
ln(sigma(h)) = ln(h^p);
ln(sigma(h)) = p*ln(h);
```
Теперь построим график с логарифмом:

![](https://github.com/timattt/Project-computational-math/blob/master/Images/GraphDiff.png)

Теперь для каждой линии p - это просто тангенс угла наклона этой линии.
И если посчитать это самое p по графику, то видно, что для первых порядков p = 1, а для вторых p = 2. Т.е. теоретическое предположени выполняется.

### Замечание
В этой задаче нужно использовать при дифференцировании симметричные формулы. Иначе в некоторых ситуациях погрешность может быть иной.   
И вообще лучше всегда во всех задачах брать все симметричное.

## Решение СЛАУ

### Метод Гаусса

Классический метод решения СЛАУ. Рассматриваем невырожденные матрица.
Про сам метод читаем [тут](https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%93%D0%B0%D1%83%D1%81%D1%81%D0%B0).

### Метод Гаусса с выбором главного элемента

Очень простая модификация предыдущего. Теперь когда обрабатываем строку выбираем такую, чтобы элемент по модулю был наибольшим.
Подроьнее [здесь](http://www.e-biblio.ru/book/bib/02_estestv_nauki/Vychislit_matematika/pr/docs/piece010.htm).

### Метод Якоби

Простейший итерационный метод

Все нужные формулы ниже:

![](https://github.com/timattt/Project-computational-math/blob/master/Images/JacobTeor.png)

D означает матрицу, у которой на главной диагонали стоят соответствующие элементы матрицы A.

Метод сходится, если 
```
||B|| < 1
```

### Метод Зейделя

Тоже итерационный метод. Возможно, немного лучше предыдущего, но все равно так себе.

Нужные соотношения ниже:

![](https://github.com/timattt/Project-computational-math/blob/master/Images/ZeidelTeor.png)

Метод сходится, если 
```
||-(L + D)^-1 * U || < 1 
```
где **L** - нижняя треугольная часть. **U** - верхняя треугольная часть. **D** - аналогично предыдущему методу. 

### Метод наискорейшего спуска

Самый потрясный итерационный метод для решения СЛАУ.

Теорминимум:

```
Дано:
Ax = f

Итерирование:
rn = Axn - f
Tn = (rn * rn) / (Arn * rn)
xn+1 = xn - Tn*rn

где Tn - число, A - входная матрица, rn - невязка
```

### Метод наименьших невязок

Аналогично предыдущему только формула для Tn другая.
```
Tn = (rn * Arn) / (Arn * Arn)
```

### Метод прогонки

Метод для решения СЛАУ, матрица которых заполнена не нулями только на трех диагоналях.

![](https://github.com/timattt/Project-computational-math/blob/master/Images/sweep.png)

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

## Интерполяции

### Метод лагранжа

Пусть у нас есть набор точек Xi и Fi на плоскости. Тогда через них можно построить полином и причем только один.

![image](https://user-images.githubusercontent.com/25401699/142418585-b8c0c74f-1e59-448f-96c0-5981f8ef7635.png)

Рассмотрим пример:

![image](https://user-images.githubusercontent.com/25401699/142419809-e843eed9-ff60-4990-9774-6138e4d86aa5.png)

### Метод Ньютона

Все аналогично предыдущему методу, но теперь выражение позволяет добавлять новые точки не пересчитывая предыдущие.

![image](https://user-images.githubusercontent.com/25401699/142418787-9f2dbc2c-ab68-431b-9126-1742c59b60b7.png)

Где разделенная разность задается выражением:

![image](https://user-images.githubusercontent.com/25401699/142418922-51d5804f-dcf1-4722-a31b-53b3a0813235.png)

Очевидно, что построенный полином будет таким же, как и в прошлом методе.
Но здесь мы можем легко добавлять еще точки. Добавим еще одну:

![image](https://user-images.githubusercontent.com/25401699/142420218-1a3902f9-f211-4f21-8c45-d3c920241ba9.png)

### Кубические сплайны

Имеем n точек, хотим построить n+1 кривых третьего порядка, чтобы их края проходили через наши точки, а производные 1 и 2 порядка не имели разрыва.

Хотим построить такие полиномы:

![image](https://user-images.githubusercontent.com/25401699/142419271-643d5a9e-3673-4de4-a7cc-ef67bd86b0cc.png)

Учитывая граничные условия:

![image](https://user-images.githubusercontent.com/25401699/142419316-4a6d5d04-f380-4c43-afd0-4d2229883140.png)

Из них получим, что

![image](https://user-images.githubusercontent.com/25401699/142419370-057bbd54-f7df-43ac-adea-56be310f3fd6.png)

Если учтем, что

![image](https://user-images.githubusercontent.com/25401699/142419425-7d7282f8-a340-4463-b265-21cd82752557.png)

то вычисление коэффициента c сведется к задаче о решении СЛАУ в трех-диагональной матрице. Решим ее методом прогонки.

Получим вот что-то такое:

![image](https://user-images.githubusercontent.com/25401699/142420356-fd41a804-0101-4e55-ad2a-09759bededac.png)

### Аппроксимация Паде

Мы хотим представить заданную функцию в виде рациональной функции.

**Алгоритм:**

![image](https://user-images.githubusercontent.com/25401699/143291967-21aec75a-a554-4f39-9d29-13c673c5fcac.png)

![image](https://user-images.githubusercontent.com/25401699/143292001-94f194d0-b88b-4211-9bfd-258667c2c7fd.png)

## Интегрирование

### Интегрирование по монтекарло

