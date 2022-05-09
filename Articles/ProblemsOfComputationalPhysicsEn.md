# Selected problems of computational physics

Some processes are constantly taking place around us. Physics gives us a series of equations that allow us to describe and predict the evolution of such phenomena.
By the way, in a wonderful book by the Strugatsky brothers, science is mentioned - the equations of mathematical magic - to describe magical processes. 
However, all these partial differential equations. And solving them is quite non-trivial. For an analytical solution, we can use methods from the Equations of Mathematical Physics (or magic).

However, our computing power is finite. Moreover, there are tasks that do not have an analytical solution at all. Therefore, computational physics comes to our aid.
By approximating smooth and continuous functions with their discrete analogues, we can teach a computer to solve almost any complex problems with good accuracy.

## Equation of thermal conductivity

The simplest example of a physical phenomenon is the propagation of heat. Let's set the task formally. Let's say we have a metal plate. It is cold at T = 0.
And let one edge of it be heated to a temperature T = 100. The question is how the templo will spread and what will be the temperature in different parts of the plate.

Let's break the plate into a grid. And we will look for the temperature values in each node of the grid. 

The equation of thermal conductivity in general looks like this:

```
Ut = Uxx + Uyy
```

The conditions for the initial temperature and for the value of the temperature at the boundary turn into the following equations:

```
U(t = 0, x, y) = U0(x, y)
U(x = 0, y) = 100
U(x = L, y) = U(x, y = 0) = U(x, y = L) = 0
```

We suddenly got a Cauchy problem. Let's now move from continuous derivatives to their numerical approximations.
For example:

```
      Ui+1 - Ui
Ut = -----------
         tau
```

There are grid nodes in the numerator here. And the denominator is the time step of the grid. Similarly, but a little more complicated, we approximate other derivatives.
Since we have 3 dimensions: x, y, t - hence the grid will be three-dimensional. Moreover, the values on all faces except one are already known from the initial conditions.

But then the difference entry in the grid looks like this:

![image](https://user-images.githubusercontent.com/25401699/166493965-8dde301d-550d-4712-8ed9-8000a8f278da.png)

And the points whose value we already know from the initial conditions are marked here with color:

![image](https://user-images.githubusercontent.com/25401699/166497773-d9f67350-fb84-4001-80ee-749b75b467a6.png)

So we have a formula to find the value at unknown points.

It remains only to ask the computer to count.

As a result, we get a qualitative dependence on time and the temperature value at all points of the rod. It can be beautifully drawn. 

![](https://user-images.githubusercontent.com/25401699/164061165-eef819fd-8f03-4a26-885d-fabdea7b28a0.gif)

This was the simplest example of using numerical methods to predict the trajectories of heat propagation. 

## String oscillation equation

Let's have a guitar now. And she has a string. And we took and pulled the string and then let go. The question is - how will the fluctuations occur?
The wave equation can answer this question. In general , it looks like this:

```
Utt = Uxx
```

Again, let's think about what we know. The ends of our string are attached to the guitar - they're not going anywhere. Moreover, we know how we pulled the string at the initial moment of time. And finally, we know that at the initial moment of time the string was at rest - after all, we held it in our hands.

These remarks can be written mathematically.

```
U(t = 0, x) = sin(x*pi/L)
U(t, x = 0) = U(t, x = L) = 0
U't(t = 0, x) = 0
```

Now let's introduce a grid. We have two variables: x and t. Therefore, the grid will be two-dimensional.
Let's mark the points on it that we already know.

![image](https://user-images.githubusercontent.com/25401699/166503975-be3883d6-b8d3-4c93-b74d-5c827777b423.png)

We also approximate numerically the derivatives. Since there is a second order, therefore, three points will already be needed, so the scheme is a cross, it is indicated in red.

Now, it is obvious that the scheme and initial conditions can be used to calculate the values at all other points.

Let's draw a final picture of what happened.

![Wave_equation_1D_fixed_endpoints](https://user-images.githubusercontent.com/25401699/166504463-e3adf09b-1bcb-48b6-ad97-3bb72d9e20df.gif)

As you can see, the wave propagates along the string, and the ends are fixed. Everything is as we ordered.
