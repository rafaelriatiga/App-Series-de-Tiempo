# ============================================================
# Paquetes del proyecto
# Se carga desde el análisis y desde la app con:
#   source("R/paquetes.R")
# ============================================================

# --- Metodología Box-Jenkins  ---
library(astsa)       # sarima(): estimación y diagnóstico
library(forecast)    # auto.arima(), ets(), forecast(), accuracy()
library(tseries)     # pruebas de raíz unitaria (adf.test)
library(nortest)     # lillie.test(): normalidad (Kolmogorov-Smirnov)
library(fBasics)     # jarqueberaTest(): normalidad
library(car)         # powerTransform() y qqPlot()
library(TSA)         # McLeod.Li.test(): efectos ARCH
library(pastecs)     # buysbal(): tabla Buys-Ballot
library(MLmetrics)   # MAPE()
library(lmtest)      # coeftest(): significancia de coeficientes

# --- App ---
library(shiny)
library(bslib)       # temas visuales para Shiny
library(plotly)      # gráficos interactivos