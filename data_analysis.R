###--------------------------------------------------------------------------###
# Código para el análisis de datos descargados desde cbioportal
###--------------------------------------------------------------------------###

### Intalación de librerias necesarias para el análisis:
if(!require("cbioportalR")){install.packages("cbioportalR")}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("gtsummary")){install.packages("gtsummary")}
if(!require("survival")){install.packages("survival")}
if(!require("survminer")){install.packages("survminer")}
if(!require("ggsurvfit")){install.packages("ggsurvfit")}
if(!require("corrplot")){install.packages("corrplot")}
if(!require("caret")){install.packages("caret")}
if(!require("cluster")){install.packages("cluster")}
if(!require("factoextra")){install.packages("factoextra")}

# Carga de datos crusdos almacenados previamente en archivo RData
load("data.RData")

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

###--------------------------------------------------------------------------###
# Analisis descriptivo de las variables
###--------------------------------------------------------------------------###
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
clinical_brca_wide$PFS_STATUS <- factor(clinical_brca_wide$PFS_STATUS, 
                                        levels = c("0:CENSORED", "1:PROGRESSION"), 
                                        labels = c("0", "1"))

# Verificar los nuevos niveles
levels(clinical_brca_wide$PFS_STATUS)

#ALTERNATIVA paquete  "ggsurvfit" --> mejor este se puede editar con ggplot2
pfs <- survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ 1, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1, color = "#607B8B") +  # Cambiar el color de la línea
  add_confidence_interval(fill = "#B0E2FF") +  # Cambiar el color del sombreado que representa el CI
  add_risktable() +  # Cambiar el estilo de la tabla de riesgo
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) 

pfs


# Zoom de los 10 primeros años (hasta 120 meses)
pfs_120m <- survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ 1, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1, color = "#607B8B") +  # Cambiar el color de la línea
  add_confidence_interval(fill = "#B0E2FF") +  # Cambiar el color del sombreado que representa el CI
  add_risktable(limits = c(0, 120)) +  # Añadir la tabla de riesgo con límites ajustados
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x = "Time (months)",
    title = "PFS"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) +  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x
  
pfs_120m




### PFS por subtipos
pfs_subtype <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2)  # Añadir marcas de censura

pfs_subtype


### PFS por subtipos  Zoom de los 10 primeros años (hasta 120 meses)
pfs_subtype_120m <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) +  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x



pfs_subtype_120m


### PFS por subtipos  Zoom de los 10 primeros años (hasta 3 años = 36 meses)
pfs_subtype_36m <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) +  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 36))  # Ajustar los límites del eje x



pfs_subtype_36m



### Análisis exploratorios: Correlación de variables numéricas
clinical_brca_numeric <- clinical_brca_wide %>% 
  select_if(is.numeric)

# lower case:
names(clinical_brca_numeric) <- c(
  "age",
  "buffa_hypoxia_score",
  "days_last_followup",
  "days_birth",
  "days_pathologic_diag",
  "dfs_month",
  "dss_month",
  "os_month",
  "pfs_month", 
  "ragnum_hypoxia_score",
  "winter_hypoxia_score"
  )

# Visualización general con un correlograma:
pairs(clinical_brca_numeric)


# Variable "days_pathologic_diag" gneera problema al calcular correlación. Eliminar:
plot(clinical_brca_numeric$days_pathologic_diag)
clinical_brca_numeric_2 <- clinical_brca_numeric
clinical_brca_numeric_2$days_pathologic_diag <- NULL
pairs(clinical_brca_numeric_2, pch=20)

# Sumario estadístico de las variables
summary(clinical_brca_numeric_2)

# Gráficos de distribuciones
par(mfrow=c(2,5))
hist(clinical_brca_numeric_2$age, main="Age")
hist(clinical_brca_numeric_2$buffa_hypoxia_score, main="Buffa hypoxia score")
hist(clinical_brca_numeric_2$days_last_followup, main="Days last followup")
hist(clinical_brca_numeric_2$days_birth, main="Days birth")
hist(clinical_brca_numeric_2$dfs_month, main="Dfs month")
hist(clinical_brca_numeric_2$dss_month, main="Dss month")
hist(clinical_brca_numeric_2$os_month, main="Os month")
hist(clinical_brca_numeric_2$pfs_month, mani="Pfs month")
hist(clinical_brca_numeric_2$ragnum_hypoxia_score, main="Ragnum hypoxia score")
hist(clinical_brca_numeric_2$winter_hypoxia_score, main="Winter hypoxia score")

# Omitir celdas con valores NA
clinical_brca_numeric_3 <- na.omit(clinical_brca_numeric_2)
summary(clinical_brca_numeric_3)
dim(clinical_brca_numeric_2)
dim(clinical_brca_numeric_3)

# Gráficos de distribuciones sin valores NAs:
par(mfrow=c(2,5))
hist(clinical_brca_numeric_3$age, main="Age")
hist(clinical_brca_numeric_3$buffa_hypoxia_score, main="Buffa hypoxia score")
hist(clinical_brca_numeric_3$days_last_followup, main="Days last followup")
hist(clinical_brca_numeric_3$days_birth, main="Days birth")
hist(clinical_brca_numeric_3$dfs_month, main="Dfs month")
hist(clinical_brca_numeric_3$dss_month, main="Dss month")
hist(clinical_brca_numeric_3$os_month, main="Os month")
hist(clinical_brca_numeric_3$pfs_month, mani="Pfs month")
hist(clinical_brca_numeric_3$ragnum_hypoxia_score, main="Ragnum hypoxia score")
hist(clinical_brca_numeric_3$winter_hypoxia_score, main="Winter hypoxia score")
# No se observa una gran cambio en las distribuciones al quitar las filas con NA

# Análisis de correlación entre todas las columnas
clinical_corr <- cor(clinical_brca_numeric_3)
corrplot(clinical_corr, type = "upper")

# Normalización de los datos
data_matrix <- as.matrix(clinical_brca_numeric_3)
normalized_matrix <- apply(data_matrix, 2, 
                           function(x) (x - min(x)) / (max(x) - min(x)))

cor_normalized <- cor(normalized_matrix)
corrplot(cor_normalized, type = "upper")

### Segcionamineto de pacientes en vase a tabla de datos clínicos
distM <- dist(normalized_matrix)
clus <- hclust(distM)
print(clus)
plot(clus)

# PCA
pca_result <- prcomp(normalized_matrix, center = TRUE, scale. = TRUE)
summary(pca_result)
plot(pca_result, type = "l",
     main="Varianza explicada por componentes principales")
library(ggfortify)
pca_plot <- autoplot(pca_result, cata = normalized_matrix)
pca_plot

biplot(pca_result)

# Revisar el otro método para hacer PCAs y sus plots

# K-means
kmeans_model <- kmeans(normalized_matrix, centers = 3)
summary(kmeans_model)
fviz_cluster(kmeans_model, data = normalized_matrix, 
             geom = c("point", "text"), main = paste("Clusters de pacientes con k-means"))

###--------------------------------------------------------------------------###
# Anslisis de datos genéticos
###--------------------------------------------------------------------------###
mutations <- genetics$mutation
cna <- genetics$cna
structural_variant <- genetics$structural_variant

# Analisis mutations 
mutations$mutationType <- as.factor(mutations$mutationType)

# Calcular el número total de mutaciones, el número de muestras únicas con al menos una mutación,
# y el porcentaje de estas muestras por gen
summary_data <- mutations %>%
  group_by(hugoGeneSymbol) %>%
  summarise(
    total_mutations = n(),  # Número total de mutaciones en este gen
    unique_samples = n_distinct(sampleId),  # Número de muestras únicas con al menos una mutación en este gen
    total_samples = n_distinct(mutations$sampleId)  # Número total de muestras en todo el dataset
  ) %>%
  mutate(
    percent_samples_with_mutation = (unique_samples / total_samples) * 100  # Porcentaje de muestras con al menos una mutación en este gen
  ) %>%
  arrange(desc(percent_samples_with_mutation)) %>%
  mutate(percent_samples_with_mutation = round(percent_samples_with_mutation, 2))



library(gt)

table_mutations <- summary_data %>%
  gt() %>%
  cols_hide(columns = "total_samples")

table_mutations
