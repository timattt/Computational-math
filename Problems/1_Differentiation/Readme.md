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

![image](https://user-images.githubusercontent.com/25401699/210582524-5ba8352d-f283-4344-b419-27970a46b06d.png)

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

![image](https://user-images.githubusercontent.com/25401699/210582639-b9fc6b0e-10d2-4879-b6db-33d7148e40fa.png)

По оси X - мелкость разбиения. По оси Y - погрешность.
Можно заметить, что у разной точности разный изгиб линии.
Чтобы найти p, прологарифмируем обе стороны:
```
sigma(h) = h ^ p;
ln(sigma(h)) = ln(h^p);
ln(sigma(h)) = p*ln(h);
```
Теперь построим график с логарифмом:

![image](https://user-images.githubusercontent.com/25401699/210582716-ba44fd76-b4ae-47ce-b2e0-4781e287a2e8.png)


Теперь для каждой линии p - это просто тангенс угла наклона этой линии.
И если посчитать это самое p по графику, то видно, что для первых порядков p = 1, а для вторых p = 2. Т.е. теоретическое предположени выполняется.

### Замечание
В этой задаче нужно использовать при дифференцировании симметричные формулы. Иначе в некоторых ситуациях погрешность может быть иной.   
И вообще лучше всегда во всех задачах брать все симметричное.
