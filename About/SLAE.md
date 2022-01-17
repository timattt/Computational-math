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