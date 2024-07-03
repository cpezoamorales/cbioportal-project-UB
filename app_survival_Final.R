library(shiny)
library(survminer)
library(survival)
library(dplyr)
library(ggplot2)
library(plotly)
library(ggsurvfit)
library(shinyWidgets)

load("data2.RData")

# Cambiar los niveles de la variable STATUS
clinical_brca_wide$PFS_STATUS <- factor(clinical_brca_wide$PFS_STATUS, 
                                        levels = c("0:CENSORED", "1:PROGRESSION"), 
                                        labels = c("0", "1"))

clinical_brca_wide$OS_STATUS <- factor(clinical_brca_wide$OS_STATUS, 
                                       levels = c("0:LIVING", "1:DECEASED"), 
                                       labels = c("0", "1"))

clinical_brca_wide$DSS_STATUS <- factor(clinical_brca_wide$DSS_STATUS, 
                                        levels = c("0:ALIVE OR DEAD TUMOR FREE", "1:DEAD WITH TUMOR"), 
                                        labels = c("0", "1"))

# Mapeo de variables a nombres amigables
variable_names <- list(
  SUBTYPE = "Subtipo",
  group_stage = "Estadío",
  age_group = "Edad"
)

# Código para crear la gráfica
plot_data <- function(time_range, variable) {
  
  # Crear gráficas para PFS, OS y DSS
  pfs_formula <- as.formula(paste("Surv(PFS_MONTHS, PFS_STATUS == '1') ~", variable))
  os_formula <- as.formula(paste("Surv(OS_MONTHS, OS_STATUS == '1') ~", variable))
  dss_formula <- as.formula(paste("Surv(DSS_MONTHS, DSS_STATUS == '1') ~", variable))
  
  pfs_plot <- survfit2(pfs_formula, data = clinical_brca_wide) %>%
    ggsurvfit(linewidth = 1) + 
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
      y = "Progression free survival",
      x = "Time (months)",
      title = paste("PFS by", variable)
    ) +
    add_censor_mark(color = "#607B8B", shape = 3, size = 1) +
    scale_x_continuous(limits = c(time_range[1], time_range[2]), expand = c(0, 0), breaks = seq(0, time_range[2], by = 24)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(
      axis.title.x = element_text(margin = margin(t = 10)), 
      axis.title.y = element_text(margin = margin(r = 10))   
    )
  
  os_plot <- survfit2(os_formula, data = clinical_brca_wide) %>%
    ggsurvfit(linewidth = 1) + 
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
      y = "Overall survival",
      x = "Time (months)",
      title = paste("OS by", variable)
    ) +
    add_censor_mark(color = "#607B8B", shape = 3, size = 1) +
    scale_x_continuous(limits = c(time_range[1], time_range[2]), expand = c(0, 0), breaks = seq(0, time_range[2], by = 24)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(
      axis.title.x = element_text(margin = margin(t = 10)), 
      axis.title.y = element_text(margin = margin(r = 10))   
    )
  
  dss_plot <- survfit2(dss_formula, data = clinical_brca_wide) %>%
    ggsurvfit(linewidth = 1) + 
    add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
    scale_ggsurvfit() +
    labs(
      y = "Disease Specific survival",
      x = "Time (months)",
      title = paste("DSS by", variable)
    ) +
    add_censor_mark(color = "#607B8B", shape = 3, size = 1) +
    scale_x_continuous(limits = c(time_range[1], time_range[2]), expand = c(0, 0), breaks = seq(0, time_range[2], by = 24)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(
      axis.title.x = element_text(margin = margin(t = 10)), 
      axis.title.y = element_text(margin = margin(r = 10))   
    )
  
  list(pfs_plot = pfs_plot, os_plot = os_plot, dss_plot = dss_plot)
}

# Definición de la interfaz de usuario
ui <- fluidPage(
  tags$head(
    tags$style(HTML('
      
      /* Estilo para sliderInput */
      .js-irs-0 .irs-bar {
        background-color: #20B2AA !important;
        box-shadow: none !important;
      }
      
      .js-irs-0 .irs-single, .js-irs-0 .irs-from, .js-irs-0 .irs-to {
        background-color: #20B2AA !important;
        color: white !important;
        font-weight: bold !important;
      }
      
      /* Ajustar el color y estilo de los números dentro del slider */
      .js-irs-0 .irs-grid-text {
        color: black !important;
        font-weight: bold !important;
      }
      
    '))
  ),
  fluidRow(
    column(12,
           h1("Visualización de curvas supervivencia en cáncer de mama", align = "center"),
           h4(style = "text-align: center; font-size: 16px;",
              "Base de datos pública, accesible en: ",
              tags$a(href = "https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018",
                     "https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018",
                     style = "color: #20B2AA; text-decoration: none;",
                     target = "_blank")
           )
    )
  ),
  fluidRow(
    column(8, offset = 2,
           sliderInput("time_range", "Rango de Tiempo (meses):",
                       min = 0, max = 240, value = c(0, 240), step = 12,
                       width = "100%", ticks = TRUE, animate = FALSE)
    ),
    style = "margin-top: 20px;"  
  ),
  fluidRow(
    column(4, offset = 4,
           selectInput("variable", "Variable de Estratificación:",
                       choices = c("Subtipo" = "SUBTYPE", "Estadío" = "group_stage", "Edad" = "age_group"))
    ),
    style = "margin-top: 20px;" 
  ),
  fluidRow(
    column(4, plotlyOutput("pfs_plot", height = "400px")),
    column(4, plotlyOutput("os_plot", height = "400px")),
    column(4, plotlyOutput("dss_plot", height = "400px")),
    style = "margin-top: 20px;"  
  ))

# Definición de la lógica del servidor
server <- function(input, output) {
  # Utilizar reactive para almacenar los gráficos generados
  plots <- reactive({
    plot_data(input$time_range, input$variable)
  })
  
  output$pfs_plot <- renderPlotly({
    ggplotly(plots()$pfs_plot, tooltip = "text") %>%
      layout(showlegend = FALSE)  # Desactiva la leyenda para este gráfico
  })
  
  output$os_plot <- renderPlotly({
    ggplotly(plots()$os_plot, tooltip = "text") %>%
      layout(showlegend = FALSE)  # Desactiva la leyenda para este gráfico
  })
  
  output$dss_plot <- renderPlotly({
    gg <- ggplotly(plots()$dss_plot, tooltip = "text")  # Convertir a plotly
    # Configurar la leyenda basada en la variable seleccionada
    gg <- gg %>% layout(showlegend = !is.null(input$variable))
    gg
  })
}
# Ejecución de la aplicación Shiny
shinyApp(ui = ui, server = server)