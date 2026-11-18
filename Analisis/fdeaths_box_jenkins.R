# ============================================================
# Aplicación SARIMA - Serie fdeaths
# Metodología Box-Jenkins (basada en el notebook 5 del curso)
# Autor: Rafael Alexander Riatiga Durán
# ============================================================

source("R/paquetes.R")

# ------------------------------------------------------------
# 1. DATOS Y ANÁLISIS DESCRIPTIVO
# ------------------------------------------------------------
# fdeaths: muertes mensuales por enfermedades pulmonares
# (bronquitis, enfisema, asma) en mujeres, Reino Unido, 1974-1979.
# Viene incluida en R (paquete datasets), no hay que cargar archivo.

y <- fdeaths

# Estructura de la serie
start(y); end(y); frequency(y); length(y)
summary(y)

# Gráfico de la serie
ts.plot(y, ylab = "Muertes", lwd = 2,
        main = "Muertes mensuales por enfermedad pulmonar - Mujeres (UK)")
grid(col = "gray")

# Gráfico con líneas verticales cada año (para ver el ciclo de 12 meses)
ts.plot(y, ylab = "Muertes", lwd = 2, main = "fdeaths con divisiones anuales")
abline(v = 1974:1980, col = "red", lty = 2)

# Tabla Buys-Ballot: filas = años, columnas = meses
tb_bb <- buysbal(y, frequency = 12)
tb_bb

# Gráficos estacionales
monthplot(y, main = "Perfil por mes (línea = promedio del mes)")
seasonplot(y, main = "Perfil estacional por año", col = 1:6, year.labels = TRUE)
boxplot(y ~ cycle(y), names = month.abb,
        main = "Distribución por mes", xlab = "Mes", ylab = "Muertes")


# Promedio anual: si hubiera tendencia, estos valores subirían o bajarían
aggregate(y, FUN = mean)


# CONCLUSIÓN 1: No hay tendencia marcada (promedios anuales entre 525 y 589).
# Estacionalidad clara de periodo 12: picos dic-mar, valles jul-sep.
# Mayor variabilidad en meses de invierno -> revisar Box-Cox.
# Dato atípico: feb-1976 (~1141).

# ------------------------------------------------------------
# 2. TRANSFORMACIÓN DE BOX-COX
# ------------------------------------------------------------
# Busca un lambda tal que y^lambda tenga varianza estable.
#   lambda = 1  -> no transformar
#   lambda = 0  -> logaritmo
#   lambda = 0.5 -> raíz cuadrada

# Alternativa 1 (car): estima lambda y hace pruebas de hipótesis
bc <- powerTransform(y)
summary(bc)

# Alternativa 2 (MASS): gráfico de verosimilitud según lambda
MASS::boxcox(y ~ 1, lambda = seq(-2, 2, 1/10), plotit = TRUE,
             xlab = expression(lambda), ylab = "Log-verosimilitud")


# Lambda estimado = -0.78, IC95% (-1.66, 0.10).
# H0: lambda = 1 rechazada (p < 0.001); H0: lambda = 0 no rechazada (p = 0.078).
# Se aplica transformación logarítmica.

yt <- log(y)

par(mfrow = c(2, 1))
plot(y,  main = "Serie original", lwd = 2); grid()
plot(yt, main = "Logaritmo de la serie", lwd = 2); grid()
par(mfrow = c(1, 1))

# ------------------------------------------------------------
# 3. IDENTIFICACIÓN: ACF Y PACF
# ------------------------------------------------------------

# Serie transformada
par(mfrow = c(2, 1))
acf(yt,  lag.max = 36, main = "ACF log(fdeaths)")
pacf(yt, lag.max = 36, main = "PACF log(fdeaths)")

# Diferencia ordinaria (1-B)
acf(diff(yt),  lag.max = 36, main = "ACF (1-B) log(fdeaths)")
pacf(diff(yt), lag.max = 36, main = "PACF (1-B) log(fdeaths)")

# Diferencia estacional (1-B^12)
acf(diff(yt, lag = 12),  lag.max = 36, main = "ACF (1-B^12) log(fdeaths)")
pacf(diff(yt, lag = 12), lag.max = 36, main = "PACF (1-B^12) log(fdeaths)")

# Ambas diferencias (1-B)(1-B^12)
acf(diff(diff(yt), lag = 12),  lag.max = 36, main = "ACF (1-B)(1-B^12) log(fdeaths)")
pacf(diff(diff(yt), lag = 12), lag.max = 36, main = "PACF (1-B)(1-B^12) log(fdeaths)")
par(mfrow = c(1, 1))

# Número de diferencias sugerido por pruebas formales
ndiffs(yt)    # ordinarias (prueba KPSS)
nsdiffs(yt)   # estacionales

# Serie log: ACF con picos en 12, 24, 36 de decaimiento lento -> estacionalidad no estacionaria.
# (1-B): no elimina el patrón estacional.
# (1-B^12): ACF dentro de bandas salvo pico negativo en rezago 12 -> MA estacional (Q = 1).
# (1-B)(1-B^12): pico negativo en rezago 1 -> sobrediferenciación; d = 0.
# Candidato: SARIMA(p,0,q)(0,1,1)_12 con p, q pequeños.
# ndiffs = 0, nsdiffs = 1: confirma d = 0, D = 1.

# ------------------------------------------------------------
# 4. IDENTIFICACIÓN POR CRITERIOS DE INFORMACIÓN
# ------------------------------------------------------------

d <- 0; D <- 1
ic <- data.frame()

for (p in 0:2) for (q in 0:2) for (P in 0:1) for (Q in 0:1) {
  fit <- tryCatch(
    Arima(yt, order = c(p, d, q), seasonal = c(P, D, Q), method = "ML"),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    ic <- rbind(ic, data.frame(p, q, P, Q, AIC = fit$aic, BIC = fit$bic))
  }
}

head(ic[order(ic$AIC), ], 5)
head(ic[order(ic$BIC), ], 5)

# Referencia: selección automática
auto.arima(yt)











