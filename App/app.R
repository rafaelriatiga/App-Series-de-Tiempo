#  ============================================================
# App Shiny - SARIMA y combinación de pronósticos para fdeaths
# Ejecutar desde la raíz del proyecto: shiny::runApp("App")
# ============================================================

library(shiny)
library(bslib)
library(plotly)
library(forecast)

# ------------------------------------------------------------
# Datos y objetos fijos
# ------------------------------------------------------------

y     <- fdeaths
train <- window(y, end = c(1978, 12))
test  <- window(y, start = c(1979, 1))
meses <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun",
           "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")

fecha <- function(x) as.Date(paste(floor(time(x) + 1e-6),
                                   round(cycle(x)), 1, sep = "-"))

mape <- function(real, pred) 100 * mean(abs((real - pred) / real))

# Modelos para la validación (entrenamiento 1974-1978, prueba 1979)
f_sarima <- forecast(Arima(train, order = c(0, 0, 0), seasonal = c(0, 1, 1),
                           lambda = 0, biasadj = TRUE), h = 12)
f_ets    <- forecast(ets(train), h = 12)
f_snaive <- snaive(train, h = 12)

pm <- cbind(SARIMA = f_sarima$mean, ETS = f_ets$mean, Naive = f_snaive$mean)

ecm <- c(mean(residuals(f_sarima, type = "response")^2, na.rm = TRUE),
         mean(residuals(f_ets,    type = "response")^2, na.rm = TRUE),
         mean(residuals(f_snaive)^2, na.rm = TRUE))
p_inv <- round((1 / ecm) / sum(1 / ecm), 2)

col_sarima <- "#2563eb"; col_ets <- "#16a34a"
col_naive  <- "#9ca3af"; col_comb <- "#dc2626"

# ------------------------------------------------------------
# Interfaz
# ------------------------------------------------------------

ui <- page_navbar(
  title = "fdeaths · Box-Jenkins",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,

  # 1. Serie --------------------------------------------------
  nav_panel(
    "1. Serie",
    layout_columns(
      value_box("Observaciones", "72", "Mensual, 1974-1979"),
      value_box("Promedio", round(mean(y)), "muertes por mes"),
      value_box("Máximo", max(y), "feb-1976 (atípico)"),
      fill = FALSE
    ),
    card(card_header("Muertes mensuales por enfermedad pulmonar - Mujeres (Reino Unido)"),
         plotlyOutput("g_serie", height = "380px")),
    card(card_header("Distribución por mes"), plotlyOutput("g_caja", height = "380px"))
  ),

  # 2. Identificación -----------------------------------------
  nav_panel(
    "2. Identificación",
    layout_sidebar(
      sidebar = sidebar(
        checkboxInput("log", "Transformación logarítmica (Box-Cox, λ ≈ 0)", TRUE),
        radioButtons("dif", "Diferenciación",
                     c("Ninguna" = "0_0", "(1-B)" = "1_0",
                       "(1-B¹²)" = "0_1", "(1-B)(1-B¹²)" = "1_1"),
                     selected = "0_1"),
        sliderInput("rez", "Rezagos", 12, 48, 36, step = 12),
        helpText("Seleccionada: d = 0, D = 1. La ACF de (1-B¹²) corta en el",
                 "rezago 12 -> MA estacional (Q = 1).")
      ),
      card(card_header("Serie transformada / diferenciada"), plotOutput("g_dif", height = 220)),
      layout_columns(
        card(card_header("ACF"),  plotOutput("g_acf",  height = 260)),
        card(card_header("PACF"), plotOutput("g_pacf", height = 260))
      )
    )
  ),

  # 3. Modelo y diagnóstico -----------------------------------
  nav_panel(
    "3. Modelo",
    layout_sidebar(
      sidebar = sidebar(
        h6("SARIMA(p, 0, q)(P, 1, Q)₁₂ sobre log(Y)"),
        layout_columns(
          numericInput("p", "p", 0, 0, 2), numericInput("q", "q", 0, 0, 2),
          numericInput("P", "P", 0, 0, 1), numericInput("Q", "Q", 1, 0, 1),
          col_widths = c(6, 6, 6, 6)
        ),
        helpText("Seleccionado: (0,0,0)(0,1,1)₁₂ — menor BIC.")
      ),
      layout_columns(
        value_box("AIC", textOutput("aic")),
        value_box("BIC", textOutput("bic")),
        uiOutput("vb_lb"), uiOutput("vb_arch"), uiOutput("vb_norm"),
        fill = FALSE
      ),
      layout_columns(
        card(card_header("Coeficientes"), tableOutput("coefs")),
        card(card_header("Residuales estandarizados"), plotOutput("g_res", height = 230)),
        col_widths = c(5, 7)
      ),
      layout_columns(
        card(card_header("ACF residuales"), plotOutput("g_acf_res", height = 230)),
        card(card_header("QQ-plot residuales"), plotOutput("g_qq", height = 230))
      )
    )
  ),

  # 4. Pronóstico ---------------------------------------------
  nav_panel(
    "4. Pronóstico",
    layout_sidebar(
      sidebar = sidebar(
        sliderInput("h", "Horizonte (meses)", 6, 36, 24, step = 6),
        sliderInput("nivel", "Nivel del intervalo (%)", 80, 99, 95)
      ),
      card(card_header("Pronóstico con el modelo de la pestaña 3 (corrección de sesgo)"),
           plotlyOutput("g_pron", height = "420px")),
      card(card_header("Tabla de pronósticos"), tableOutput("t_pron"))
    )
  ),

  # 5. Combinación de modelos ---------------------------------
  nav_panel(
    "5. Combinación",
    layout_sidebar(
      sidebar = sidebar(
        h6("Pesos de cada modelo"),
        sliderInput("w1", "SARIMA", 0, 1, 1/3, step = 0.01),
        sliderInput("w2", "ETS (Holt-Winters)", 0, 1, 1/3, step = 0.01),
        sliderInput("w3", "Naive estacional", 0, 1, 1/3, step = 0.01),
        helpText("Los pesos se normalizan para que sumen 1."),
        actionButton("b_igual", "Pesos iguales", class = "btn-outline-primary"),
        actionButton("b_inv", "Inverso del ECM", class = "btn-outline-primary"),
        actionButton("b_sar", "Solo SARIMA", class = "btn-outline-secondary")
      ),
      layout_columns(
        uiOutput("vb_comb"), uiOutput("vb_mejor"), uiOutput("vb_pesos"),
        fill = FALSE
      ),
      card(card_header("Validación 1979: real vs pronósticos (entrenamiento 1974-1978)"),
           plotlyOutput("g_comb", height = "420px")),
      card(card_header("MAPE en 1979"), tableOutput("t_mape"))
    )
  )
)

# ------------------------------------------------------------
# Servidor
# ------------------------------------------------------------

server <- function(input, output, session) {

  # 1. Serie --------------------------------------------------
  output$g_serie <- renderPlotly({
    plot_ly(x = fecha(y), y = as.numeric(y), type = "scatter", mode = "lines+markers",
            line = list(color = col_sarima), marker = list(size = 5),
            hovertemplate = "%{x|%b %Y}: %{y}<extra></extra>") |>
      layout(xaxis = list(title = ""), yaxis = list(title = "Muertes"))
  })

  output$g_caja <- renderPlotly({
    plot_ly(x = factor(meses[cycle(y)], levels = meses), y = as.numeric(y),
            type = "box", marker = list(color = col_sarima),
            line = list(color = col_sarima)) |>
      layout(xaxis = list(title = ""), yaxis = list(title = "Muertes"))
  })

  # 2. Identificación -----------------------------------------
  serie_id <- reactive({
    z <- if (input$log) log(y) else y
    d <- as.numeric(strsplit(input$dif, "_")[[1]])
    if (d[2] == 1) z <- diff(z, lag = 12)
    if (d[1] == 1) z <- diff(z)
    z
  })

  output$g_dif  <- renderPlot({ par(mar = c(3, 4, 1, 1)); plot(serie_id(), ylab = "", lwd = 2); grid() })
  output$g_acf  <- renderPlot({ par(mar = c(4, 4, 1, 1)); acf(serie_id(),  lag.max = input$rez, main = "") })
  output$g_pacf <- renderPlot({ par(mar = c(4, 4, 1, 1)); pacf(serie_id(), lag.max = input$rez, main = "") })

  # 3. Modelo -------------------------------------------------
  modelo <- reactive({
    req(input$p, input$q, input$P, input$Q)
    tryCatch(
      Arima(y, order = c(input$p, 0, input$q), seasonal = c(input$P, 1, input$Q),
            lambda = 0, biasadj = TRUE, method = "ML"),
      error = function(e) NULL
    )
  })

  residuos <- reactive({ req(modelo()); residuals(modelo()) })

  output$aic <- renderText({ req(modelo()); round(modelo()$aic, 2) })
  output$bic <- renderText({ req(modelo()); round(modelo()$bic, 2) })

  caja_prueba <- function(titulo, p, ok_si_mayor = TRUE) {
    ok <- if (ok_si_mayor) p > 0.05 else p < 0.05
    p_txt <- if (p < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p))
    value_box(titulo, p_txt,
              if (ok) "Supuesto cumplido" else "Supuesto no cumplido",
              theme = if (ok) "success" else "danger")
  }

  output$vb_lb <- renderUI({
    k <- length(coef(modelo()))
    caja_prueba("Ljung-Box (incorrelación)",
                Box.test(residuos(), lag = 24, type = "Ljung-Box", fitdf = k)$p.value)
  })
  output$vb_arch <- renderUI({
    caja_prueba("McLeod-Li (sin ARCH)",
                Box.test(residuos()^2, lag = 12, type = "Ljung-Box")$p.value)
  })
  output$vb_norm <- renderUI({
    caja_prueba("Shapiro-Wilk (normalidad)", shapiro.test(residuos())$p.value)
  })

  output$coefs <- renderTable({
    m <- modelo(); req(m)
    se <- sqrt(diag(m$var.coef))
    z  <- coef(m) / se
    data.frame("Parámetro" = names(coef(m)), "Estimación" = coef(m), "Error est." = se,
               "p-valor" = 2 * pnorm(-abs(z)), check.names = FALSE)
  }, digits = 4)

  output$g_res <- renderPlot({
    r <- residuos() / sqrt(modelo()$sigma2)
    par(mar = c(3, 4, 1, 1))
    plot(r, ylab = "", type = "h", col = col_sarima)
    abline(h = c(-3, 3), col = col_comb, lty = 2)
  })
  output$g_acf_res <- renderPlot({ par(mar = c(4, 4, 1, 1)); acf(residuos(), lag.max = 36, main = "") })
  output$g_qq <- renderPlot({
    par(mar = c(4, 4, 1, 1))
    qqnorm(residuos(), main = "", pch = 19, col = col_sarima)
    qqline(residuos(), col = col_comb, lwd = 2)
  })

  # 4. Pronóstico ---------------------------------------------
  pronostico <- reactive({
    req(modelo())
    forecast(modelo(), h = input$h, level = input$nivel)
  })

  output$g_pron <- renderPlotly({
    f <- pronostico()
    plot_ly() |>
      add_ribbons(x = fecha(f$mean), ymin = as.numeric(f$lower), ymax = as.numeric(f$upper),
                  fillcolor = "rgba(37,99,235,0.2)", line = list(width = 0),
                  name = paste0("IC ", input$nivel, "%")) |>
      add_lines(x = fecha(y), y = as.numeric(y), name = "Observado",
                line = list(color = "black")) |>
      add_lines(x = fecha(f$mean), y = as.numeric(f$mean), name = "Pronóstico",
                line = list(color = col_sarima, width = 3)) |>
      layout(xaxis = list(title = ""), yaxis = list(title = "Muertes"))
  })

  output$t_pron <- renderTable({
    f <- pronostico()
    data.frame(Mes = format(fecha(f$mean), "%b %Y"),
               "Pronóstico" = as.numeric(f$mean),
               Inferior = as.numeric(f$lower), Superior = as.numeric(f$upper),
               check.names = FALSE)
  }, digits = 1)

  # 5. Combinación --------------------------------------------
  fijar_pesos <- function(w) {
    updateSliderInput(session, "w1", value = w[1])
    updateSliderInput(session, "w2", value = w[2])
    updateSliderInput(session, "w3", value = w[3])
  }
  observeEvent(input$b_igual, fijar_pesos(c(1, 1, 1) / 3))
  observeEvent(input$b_inv,   fijar_pesos(p_inv))
  observeEvent(input$b_sar,   fijar_pesos(c(1, 0, 0)))

  pesos <- reactive({
    w <- c(input$w1, input$w2, input$w3)
    if (sum(w) == 0) w <- c(1, 1, 1)
    w / sum(w)
  })

  combinado <- reactive(ts(pm %*% pesos(), start = c(1979, 1), frequency = 12))

  tabla_mape <- reactive({
    data.frame(
      Modelo = c("SARIMA", "ETS (Holt-Winters)", "Naive estacional", "Combinación"),
      Peso = c(pesos(), 1),
      MAPE = c(mape(test, pm[, 1]), mape(test, pm[, 2]),
               mape(test, pm[, 3]), mape(test, combinado()))
    )
  })

  output$vb_comb <- renderUI({
    value_box("MAPE combinación", paste0(round(mape(test, combinado()), 2), " %"),
              theme = "primary")
  })
  output$vb_mejor <- renderUI({
    t <- tabla_mape()[1:3, ]
    value_box("Mejor modelo individual", t$Modelo[which.min(t$MAPE)],
              paste0("MAPE = ", round(min(t$MAPE), 2), " %"))
  })
  output$vb_pesos <- renderUI({
    w <- round(pesos(), 2)
    value_box("Pesos normalizados", paste(w, collapse = " · "), "SARIMA · ETS · Naive")
  })

  output$g_comb <- renderPlotly({
    x <- fecha(test)
    plot_ly() |>
      add_lines(x = fecha(window(y, start = c(1977, 1))),
                y = as.numeric(window(y, start = c(1977, 1))),
                name = "Observado", line = list(color = "black", width = 3)) |>
      add_lines(x = x, y = as.numeric(pm[, 1]), name = "SARIMA",
                line = list(color = col_sarima, dash = "dot")) |>
      add_lines(x = x, y = as.numeric(pm[, 2]), name = "ETS",
                line = list(color = col_ets, dash = "dot")) |>
      add_lines(x = x, y = as.numeric(pm[, 3]), name = "Naive est.",
                line = list(color = col_naive, dash = "dot")) |>
      add_lines(x = x, y = as.numeric(combinado()), name = "Combinación",
                line = list(color = col_comb, width = 4)) |>
      layout(xaxis = list(title = ""), yaxis = list(title = "Muertes"))
  })

  output$t_mape <- renderTable(tabla_mape(), digits = 2)
}

shinyApp(ui, server)
