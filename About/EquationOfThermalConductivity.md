# Уравнение теплопроводности

## Одномерный случай

Хотим решить вот такую вещь:

```
Ut = Uxx
```

С начальными условиями:

```
U(0, x) = phi(x)
U(t, 0) = psi(t)
U(t, L) = mu(t)
```

Итого на сетке начальные условия покрывают три грани.
Значит можно использовать для решения общую схему с сигмой.
Вот такую:

```
*---*---*
    |
*---*---*
```

Посколько из такой схемы выделить переменные адекватно нельзя, то просто сделаем одну жирную матрицу со всеми неизвестными значениями функции и решим ее.

Получим что-то такое:

![image](https://github.com/timattt/Project-computational-math/blob/master/Images/animation.gif)

![image](https://github.com/timattt/Project-computational-math/blob/master/Images/animation1.gif)

### Оценка погрешностей

![image](https://user-images.githubusercontent.com/25401699/164059355-3d0993d1-c770-4924-bddd-657ceee0761d.png)


## Двумерный случай

Используем вот такую схему:

![image](https://user-images.githubusercontent.com/25401699/164059197-bfabc474-1a2d-4d7d-9fd0-c00a76944c1c.png)

А также зададим начальные условия вот так:

![image](https://user-images.githubusercontent.com/25401699/164059304-b4c00351-58df-437d-b67c-75a2d97f8b31.png)

И получим что-то такое:

![1_8h4MOW_Oy7sXdidBJQ1p8g](https://user-images.githubusercontent.com/25401699/164061165-eef819fd-8f03-4a26-885d-fabdea7b28a0.gif)
