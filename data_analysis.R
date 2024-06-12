
library(cbioportalR)
library(dplyr)
library (tidyverse)
library(purrr)
library(gtsummary)
library(survival)
library(survminer)
library(ggsurvfit)



# Load data, this format is faster thar read from the cbioportal API
load("C:/Users/ljoval/OneDrive - VHIO/LAIA/CURSOS/UB_DATA_SCIENCE/TREBALL FINAL/cbioportal-project-UB/data.RData")


# Transformar la tabla en Wide con pivot_wider, crea columnas viniendo de muchas filas
clinical_brca_wide <- clinical_brca %>%
  select(patientId, clinicalAttributeId, value)  %>%
  pivot_wider(names_from = clinicalAttributeId, values_from = value)

str(clinical_brca_wide)

length(unique(clinical_brca$patientId))


#Transformar variables numericas a numericas y factores

convert_type <- function(column, type) {
  if (type == "NUMBER") {
    return(as.numeric(column))
  } else if (type == "STRING") {
    return(as.factor(column))
  } else {
    return(column)  # Por defecto, retornar la columna sin cambios
  }
}

# Convertir las columnas de clinical_brca_wide según el datatype de attr_brca
for (i in seq_len(nrow(attr_brca))) {
  column_name <- attr_brca$clinicalAttributeId[i]
  column_type <- attr_brca$datatype[i]
  
  if (column_name %in% names(clinical_brca_wide)) {
    clinical_brca_wide[[column_name]] <- convert_type(clinical_brca_wide[[column_name]], column_type)
  }
  
}


# Analisis descriptivo de las variables

# Crear resumen de tabla excluyendo las columnas uniquePatientKey, patientId y studyId
tbl <- tbl_summary(
  clinical_brca_wide %>%
    select(-c(patientId, FORM_COMPLETION_DATE, ICD_10, OTHER_PATIENT_ID )), 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels()

tbl

# Cambiar los niveles de la variable PFS_STATUS
clinical_brca_wide$PFS_STATUS <- factor(clinical_brca_wide$PFS_STATUS, levels = c("0:CENSORED", "1:PROGRESSION"), labels = c("0", "1"))

# Verificar los nuevos niveles
levels(clinical_brca_wide$PFS_STATUS)



# Crear un objeto de supervivencia
 pfs <- survfit(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide)

#PFS_STATUS == "1" es per generar un vector TRUE, FALSE i que agafi d'event el TRUE

pfs

# Crear el gráfico de supervivencia
ggsurvplot(pfs, data = clinical_brca_wide, risk.table = TRUE, legend.title="", title= "PFS", xlab="Time (months)", ylab= "Progression Free Survival (PFS)")


#ALTERNATIVA paquete  "ggsurvfit" --> mejor este se puede editar con ggplot2


#> Loading required package: ggplot2

p <- survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ 1, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1, color = "#607B8B") +  # Cambiar el color de la línea
  add_confidence_interval(fill = "#B0E2FF") +  # Cambiar el color del sombreado
  add_risktable() +  # Cambiar el estilo de la tabla de riesgo
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS"
  )

p


tbl_survfit(pfs, times=c(60, 120),label_header = "% Survival at Time {time} (95% CI)")


p2 <- survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1, color = "#607B8B") +  # Cambiar el color de la línea
  add_confidence_interval(fill = "#B0E2FF") +  # Cambiar el color del sombreado
  add_risktable() +  # Cambiar el estilo de la tabla de riesgo
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS"
  )

p2


