# App Series de Tiempo: SARIMA para fdeaths

Aplicación de la metodología Box-Jenkins a la serie `fdeaths` de R: muertes mensuales por enfermedades pulmonares (bronquitis, enfisema y asma) en mujeres del Reino Unido, 1974-1979 (72 observaciones). Incluye el análisis completo en R y una app interactiva en Shiny.

**Curso:** Series de Tiempo, Universidad de Antioquia
**Autor:** Rafael Alexander Riatiga Durán

## Estructura

```
├── Analisis/fdeaths_box_jenkins.R   # Análisis Box-Jenkins paso a paso
├── App/app.R                        # App Shiny
├── R/paquetes.R                     # Paquetes del proyecto
├── renv.lock                        # Versiones exactas de los paquetes
└── .Rprofile, renv/                 # Activación del entorno renv
```

## Cómo ejecutarlo

Requiere R 4.5 o superior.

```r
# 1. Clonar el repositorio y abrir R en la carpeta del proyecto.
# 2. Instalar los paquetes con las versiones del proyecto:
renv::restore()

# 3. Abrir la app:
shiny::runApp("App")
```

El análisis se corre línea por línea desde `Analisis/fdeaths_box_jenkins.R`.

## Metodología y resultados

| Etapa | Resultado |
|---|---|
| Descriptivo | Sin tendencia marcada; estacionalidad de periodo 12 (picos dic-mar, valles jul-sep). Atípico en feb-1976. |
| Box-Cox | λ = -0.78, IC95% (-1.66, 0.10). Se rechaza λ = 1 y no se rechaza λ = 0 → logaritmo. |
| Identificación | ACF/PACF de (1-B¹²): corte en el rezago 12. `ndiffs = 0`, `nsdiffs = 1` → d = 0, D = 1. |
| Selección | AIC y BIC sobre 36 modelos. Diferencia de AIC < 2 entre los mejores; se elige el de menor BIC. |
| Modelo | **SARIMA(0,0,0)(0,1,1)₁₂ sobre log(Y)**: (1 - B¹²) log Yₜ = (1 - 0.713 B¹²) aₜ, con sma1 significativo (p < 0.001). |
| Diagnóstico | Ljung-Box p = 0.81 (ruido blanco). McLeod-Li p = 0.82 (sin efectos ARCH). Normalidad rechazada por el atípico de feb-1976. |
| Pronóstico | 24 meses con retransformación corregida por sesgo: exp(m + s²/2). |

### Validación fuera de muestra y combinación de pronósticos

Entrenamiento 1974-1978 y prueba en 1979:

| Modelo | MAPE 1979 |
|---|---|
| SARIMA(0,0,0)(0,1,1)₁₂ | 5.07 % |
| Combinación (pesos inversos al ECM) | 5.43 % |
| ETS (Holt-Winters) | 5.56 % |
| Combinación (pesos iguales) | 5.64 % |
| Naive estacional | 7.23 % |

En esta serie la combinación no supera al SARIMA: cuando un modelo ya captura bien la estructura, mezclarlo con modelos menos precisos diluye su ventaja.

## App Shiny

1. **Serie:** gráfico interactivo y distribución por mes.
2. **Identificación:** transformación y diferenciación seleccionables, con su ACF/PACF.
3. **Modelo:** órdenes p, q, P, Q ajustables; AIC/BIC, coeficientes y pruebas de los residuales (verde = supuesto cumplido).
4. **Pronóstico:** horizonte y nivel del intervalo ajustables.
5. **Combinación:** pesos de SARIMA, ETS y naive estacional ajustables, con su MAPE en 1979.
