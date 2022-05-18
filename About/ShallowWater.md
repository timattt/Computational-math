# Уравнения мелкой воды

## Ссылки 

[Уравнения Навье-Стокса](https://ru.wikipedia.org/wiki/%D0%A3%D1%80%D0%B0%D0%B2%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F_%D0%9D%D0%B0%D0%B2%D1%8C%D0%B5_%E2%80%94_%D0%A1%D1%82%D0%BE%D0%BA%D1%81%D0%B0)

[Уравнения мелкой воды](https://ru.wikipedia.org/wiki/%D0%A3%D1%80%D0%B0%D0%B2%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F_%D0%BC%D0%B5%D0%BB%D0%BA%D0%BE%D0%B9_%D0%B2%D0%BE%D0%B4%D1%8B)

## Уравнение Хопфа
 
 Это простейшее нелинейное уравнение в частных производных из гидродинамики.
 
### Постановка задачи
 
 ![image](https://user-images.githubusercontent.com/25401699/167414202-728b93ed-815d-44ec-8113-2f57eed345c3.png)

Начальные условия в точности, как и в прошлой задаче.

### Решение

Используем метод характеристик.

![image](https://user-images.githubusercontent.com/25401699/167414406-e5a81de4-b9f0-460c-8bfd-8736ebd907aa.png)

Рассуждения в точности, как и в предыдущей задаче. Только нам нужно отдельно потребовать, чтобы U было больше нуля всегда.

![image](https://user-images.githubusercontent.com/25401699/167414527-c7595f8d-935c-4712-97d9-0c6f59bc50f6.png)

Итого, имеем сетку решения.

### Тестирование

![](https://github.com/timattt/Project-computational-math/blob/master/Images/hopf.gif)
