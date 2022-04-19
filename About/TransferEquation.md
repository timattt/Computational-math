# Уравнение в частных производных первого порядка

## Постановка задачи

```
Ut + a*Ux = 0
U(t = 0) = phi(x)
U(x = 0) = mu(t)
```

## Простой тестовый пример

```
a = 1
U(t = 0) = sin(x)
U(x = 0) = -sin(t)
U = sin(x - t)
```

## Схема простой уголок

```
*---*
|
*
```

![image](https://user-images.githubusercontent.com/25401699/160690160-c3dbc11f-93db-4111-b840-539bc73b195a.png)

## Схема неявный уголок

```
*---*
    |
    *
```

![image](https://user-images.githubusercontent.com/25401699/160690261-cfdf516a-0862-4bae-a34f-378e92f88794.png)

Как видим, теперь решение съезжает.

## Схема Лакса-Вендрофа

```
    *
    |
*---*---*
```

![image](https://user-images.githubusercontent.com/25401699/160690431-fa335ae1-08a1-451f-97fd-2781c4ccb41c.png)

## Оценка погрешностей

![image](https://user-images.githubusercontent.com/25401699/160690565-28039c77-eba3-49f3-879c-cdaa53071364.png)

| Простой уголок                  | Неявный уголок           | Лакс-Вендорф                |
|---------------------------------|--------------------------|-----------------------------|
| 1                               | 1                        | 2                           |

## Более прикольные примеры

### Колокол

![animation](https://user-images.githubusercontent.com/25401699/160690919-f9fed6e8-fefe-4de9-a62e-3b686c10d9b3.gif)

### Плато слева

![animation](https://user-images.githubusercontent.com/25401699/160691042-99f7849e-fd7a-4d2c-be77-64b3a677febb.gif)

### Плато справа

![animation](https://user-images.githubusercontent.com/25401699/160691203-7d30cf29-604d-4ffc-b761-ce49f3cfb706.gif)

### Плато посередине

![animation](https://user-images.githubusercontent.com/25401699/160691268-76b18a94-83c0-46c2-ae6c-955894086858.gif)

### Сшивка

Берем схему ЛВ. Замечаем, что она состоит из уголка и чего-то еще. Используем уголок. Но в моменты разрыва производной добавляем это самое что-то еще.

![animation](https://user-images.githubusercontent.com/25401699/160696843-d3e920f0-963c-4bd8-8fc8-24e6e35e64d1.gif)
