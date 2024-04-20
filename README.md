# Exponentes-de-Lyapunov
Código Fortran que implementa el método de Benettin para realizar el cálculo de los exponentes de Lyapunov y posteriormente seleccionando un parámetro del modelo ('a') realiza un espectro de los Exponentes de Lyapunov

```math
\begin{gathered}
    \overset{\cdot}{N}=c\cdot N(1-N)-d\cdot TN\\
    \overset{\cdot}{T}= T(1-T)-a\cdot TN-b\cdot IT\\
    \overset{\cdot}{I}=\frac{e \cdot IT}{f + T}-g\cdot IT -h\cdot I
\end{gathered} 
```
Este sistema es el mostrado en Itik, M., & Banks, S. P. (2010). Chaos in a three-dimensional cancer model. International Journal of Bifurcation and Chaos, 20(01), 71-79. https://doi.org/10.1142/S0218127410025417
![Espectro de los Exponentes Caracteristicos de Lyapunov](https://github.com/Cygnus000/Exponentes-de-Lyapunov/blob/main/lyapuniv-c.png)
![Espectro de los Exponentes Caracteristicos de Lyapunov](https://github.com/Cygnus000/Exponentes-de-Lyapunov/blob/main/lyapuniv-c(puntos).png)
![Espectro de los Exponentes Caracteristicos de Lyapunov](https://github.com/Cygnus000/Exponentes-de-Lyapunov/blob/main/lyapuniv-c(positivo).png)
