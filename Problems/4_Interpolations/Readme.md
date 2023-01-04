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
