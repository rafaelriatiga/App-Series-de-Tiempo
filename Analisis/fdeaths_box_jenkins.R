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

# ------------------------------------------------------------
# 5. ESTIMACIÓN
# ------------------------------------------------------------

m1 <- Arima(yt, order = c(0, 0, 0), seasonal = c(0, 1, 1), method = "ML")
m2 <- Arima(yt, order = c(0, 0, 1), seasonal = c(0, 1, 1), method = "ML")
m3 <- Arima(yt, order = c(0, 0, 0), seasonal = c(0, 1, 1),
            include.drift = TRUE, method = "ML")

coeftest(m1)
coeftest(m2)
coeftest(m3)

# m1: sma1 = -0.713 (p < 0.001), significativo.
# m2: ma1 = 0.250 (p = 0.094), no significativo -> no mejora a m1.
# m3: sma1 = -0.999 (frontera de invertibilidad, p = 0.10); el drift absorbe
#     la estacionalidad -> modelo inestable, se descarta.
# Modelo seleccionado: SARIMA(0,0,0)(0,1,1)_12 sobre log(fdeaths).
#   (1 - B^12) log(Y_t) = (1 - 0.713 B^12) a_t

mod <- m1
mod

# ------------------------------------------------------------
# 6. DIAGNÓSTICO DE RESIDUALES
# ------------------------------------------------------------

res <- residuals(mod)
res_est <- res / sqrt(mod$sigma2)

# Valores ajustados
ts.plot(yt, fitted(mod), col = c("black", "red"), lwd = 2,
        main = "log(fdeaths): observado (negro) y ajustado (rojo)")

# Residuales estandarizados (bandas en +-3)
plot(res_est, main = "Residuales estandarizados", ylab = "")
abline(h = c(-3, 0, 3), col = c("red", "gray", "red"), lty = 2)

# Incorrelación
par(mfrow = c(2, 1))
acf(res,  lag.max = 36, main = "ACF residuales")
pacf(res, lag.max = 36, main = "PACF residuales")
par(mfrow = c(1, 1))

Box.test(res, lag = 24, type = "Ljung-Box", fitdf = 1)

# Normalidad
shapiro.test(res)
lillie.test(res)
jarqueberaTest(res)
qqPlot(res, envelope = 0.95, main = "QQ-plot residuales")

# Heterocedasticidad condicional (efectos ARCH):
# Ljung-Box sobre residuales al cuadrado (prueba de McLeod-Li)
Box.test(res^2, lag = 12, type = "Ljung-Box")

# Ljung-Box p = 0.81: residuales incorrelacionados (ruido blanco).
# Ljung-Box sobre res^2 p = 0.82: sin efectos ARCH; no se requiere GARCH.
# Normalidad rechazada (Shapiro p < 0.001), asociada al atípico de feb-1976
# (residual estandarizado > 3). Afecta los intervalos, no los pronósticos puntuales.

# ------------------------------------------------------------
# 7. PRONÓSTICOS
# ------------------------------------------------------------

h <- 24
pred <- predict(mod, n.ahead = h)

# Retransformación con corrección de sesgo: exp(m + s^2 / 2)
pron   <- exp(pred$pred + 0.5 * pred$se^2)
li     <- exp(pred$pred - 1.96 * pred$se + 0.5 * pred$se^2)
ls     <- exp(pred$pred + 1.96 * pred$se + 0.5 * pred$se^2)

ts.plot(y, pron, li, ls, col = c("black", "blue", "red", "red"),
        lty = c(1, 1, 2, 2), lwd = 2,
        main = "Pronóstico SARIMA(0,0,0)(0,1,1)_12 - 24 meses")

pronosticos <- cbind(Pronostico = pron, Li_95 = li, Ls_95 = ls)
round(pronosticos, 1)

# ------------------------------------------------------------
# 8. VALIDACIÓN FUERA DE MUESTRA Y COMBINACIÓN DE MODELOS
# ------------------------------------------------------------
# Entrenamiento: 1974-1978 (60 meses). Prueba: 1979 (12 meses).

train <- window(y, end = c(1978, 12))
test  <- window(y, start = c(1979, 1))

# lambda = 0 -> log; biasadj = TRUE -> corrección de sesgo al retransformar
f_sarima <- forecast(Arima(train, order = c(0, 0, 0), seasonal = c(0, 1, 1),
                           lambda = 0, biasadj = TRUE), h = 12)
f_ets    <- forecast(ets(train), h = 12)
f_snaive <- snaive(train, h = 12)

# Pesos iguales
p_igual <- c(1, 1, 1) / 3

# Pesos inversos al error cuadrático medio de entrenamiento
ecm <- c(mean(residuals(f_sarima, type = "response")^2, na.rm = TRUE),
         mean(residuals(f_ets,    type = "response")^2, na.rm = TRUE),
         mean(residuals(f_snaive)^2, na.rm = TRUE))
p_inv <- (1 / ecm) / sum(1 / ecm)
round(p_inv, 3)

pm <- cbind(f_sarima$mean, f_ets$mean, f_snaive$mean)
f_comb_igual <- ts(pm %*% p_igual, start = c(1979, 1), frequency = 12)
f_comb_inv   <- ts(pm %*% p_inv,   start = c(1979, 1), frequency = 12)

mape <- function(real, pred) 100 * mean(abs((real - pred) / real))

tabla_mape <- data.frame(
  Modelo = c("SARIMA", "ETS (Holt-Winters)", "Naive estacional",
             "Combinación pesos iguales", "Combinación inverso ECM"),
  MAPE = round(c(mape(test, f_sarima$mean), mape(test, f_ets$mean),
                 mape(test, f_snaive$mean), mape(test, f_comb_igual),
                 mape(test, f_comb_inv)), 2)
)
tabla_mape

ts.plot(test, f_sarima$mean, f_ets$mean, f_snaive$mean, f_comb_inv,
        col = c("black", "blue", "darkgreen", "gray50", "red"),
        lwd = c(3, 2, 2, 2, 3), lty = c(1, 1, 1, 2, 1),
        main = "1979: real vs pronósticos")
legend("topright", bty = "n", lwd = 2,
       col = c("black", "blue", "darkgreen", "gray50", "red"),
       legend = c("Real", "SARIMA", "ETS", "Naive est.", "Combinación"))










