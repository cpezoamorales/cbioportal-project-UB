###--------------------------------------------------------------------------###
# Código para el análisis de datos descargados desde cbioportal
###--------------------------------------------------------------------------###

### Intalación de librerias necesarias para el análisis:
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("gtsummary")){install.packages("gtsummary")}
if(!require("survival")){install.packages("survival")}
if(!require("survminer")){install.packages("survminer")}
if(!require("ggsurvfit")){install.packages("ggsurvfit")}
if(!require("corrplot")){install.packages("corrplot")}
if(!require("caret")){install.packages("caret")}
if(!require("cluster")){install.packages("cluster")}
if(!require("factoextra")){install.packages("factoextra")}

# Carga de datos crudos almacenados previamente en archivo RData
load("data.RData")

# Transformar la tabla en Wide con pivot_wider, crea columnas viniendo de muchas filas
clinical_brca_wide <- clinical_brca %>%
  select(patientId, clinicalAttributeId, value)  %>%
  pivot_wider(names_from = clinicalAttributeId, values_from = value)

str(clinical_brca_wide)
dim(clinical_brca_wide)
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

str(clinical_brca_wide)

###-------------------------------------------------------------------------###
# Análisis descriptivo de las variables
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

###--------------------------------------------------------------------------###
### Cálculo de sobrevida de las pacientes
###--------------------------------------------------------------------------###
# Verificar los nuevos niveles
levels(clinical_brca_wide$PFS_STATUS)

p <- survival::survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ 1, data = clinical_brca_wide) %>%
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

### PFS por subtipos
pfs_subtype <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by subtype"
  )

pfs_subtype

###--------------------------------------------------------------------------###
### Análisis de correlación entre variables numéricas
###--------------------------------------------------------------------------###

# Creación de una nueva tabla que contiene solo los datos numéricos:
clinical_brca_numeric <- clinical_brca_wide %>% 
  select_if(is.numeric)

clinical_selection <- clinical_brca_wide %>% 
  select(AGE, BUFFA_HYPOXIA_SCORE, RAGNUM_HYPOXIA_SCORE, WINTER_HYPOXIA_SCORE, GENETIC_ANCESTRY_LABEL)

# Transformación de los nombres de las variables:
names(clinical_selection) <- c(
  "age",
  "buffa_hypoxia_score",
  "ragnum_hypoxia_score",
  "winter_hypoxia_score",
  "genetic_ancestry"
  )

dim(clinical_selection)
clinical_selection <- na.omit(clinical_selection)
dim(clinical_selection)

# Crear valores dummy para variabel ancestría genetica
library(fastDummies)
df <- dummy_cols(clinical_selection, select_columns = "genetic_ancestry")
clinical_selection_2 <- df[, -5]


# Visualización general con un correlograma:
pairs(clinical_brca_numeric)

# Variable "days_pathologic_diag" gneera problema al calcular correlación. Eliminar:
plot(clinical_brca_numeric$days_pathologic_diag)
clinical_brca_numeric_2 <- clinical_brca_numeric
clinical_brca_numeric_2$days_pathologic_diag <- NULL
pairs(clinical_brca_numeric_2, pch=20)

# Sumario estadístico de las variables
summary(clinical_brca_numeric_2)

# Selección de variables de interés
clinical_brca_numeric_2 <- clinical_brca_numeric_2 %>% 
  select(age, buffa_hypoxia_score, ragnum_hypoxia_score, winter_hypoxia_score)


# Gráficos de distribuciones
par(mfrow=c(1,4))
hist(clinical_brca_numeric_2$age, main="Age")
hist(clinical_brca_numeric_2$buffa_hypoxia_score, main="Buffa hypoxia score")
hist(clinical_brca_numeric_2$ragnum_hypoxia_score, main="Ragnum hypoxia score")
hist(clinical_brca_numeric_2$winter_hypoxia_score, main="Winter hypoxia score")

# Omitir celdas con valores NA
clinical_brca_numeric_3 <- na.omit(clinical_brca_numeric_2)
summary(clinical_brca_numeric_3)
dim(clinical_brca_numeric_2)
dim(clinical_brca_numeric_3)

# Gráficos de distribuciones sin valores NAs:
# Gráficos de distribuciones
par(mfrow=c(2,4))
hist(clinical_brca_numeric_2$age, main="Age")
hist(clinical_brca_numeric_2$buffa_hypoxia_score, main="Buffa hypoxia score")
hist(clinical_brca_numeric_2$ragnum_hypoxia_score, main="Ragnum hypoxia score")
hist(clinical_brca_numeric_2$winter_hypoxia_score, main="Winter hypoxia score")
hist(clinical_brca_numeric_3$age, main="Age (2)")
hist(clinical_brca_numeric_3$buffa_hypoxia_score, main="Buffa hypoxia score (2)")
hist(clinical_brca_numeric_3$ragnum_hypoxia_score, main="Ragnum hypoxia score (2)")
hist(clinical_brca_numeric_3$winter_hypoxia_score, main="Winter hypoxia score (2)")
# No se observa una gran cambio en las distribuciones al quitar las filas con NA

# Análisis de correlación entre todas las columnas
clinical_corr <- cor(clinical_brca_numeric_3)
par(mfrow=c(1,1))
corrplot(clinical_corr, type = "upper")

# Normalización de los datos
data_matrix <- as.matrix(clinical_brca_numeric_3)
normalized_matrix <- apply(data_matrix, 2, 
                           function(x) (x - min(x)) / (max(x) - min(x)))

data_matrix <- as.matrix(clinical_selection_2)
normalized_matrix <- apply(data_matrix, 2, 
                           function(x) (x - min(x)) / (max(x) - min(x)))


cor_normalized <- cor(normalized_matrix)
corrplot(cor_normalized, type = "upper")

### Clusters de pacientes en vase a tabla de datos clínicos
# Agregar ancestría como variable dummy

distM <- dist(normalized_matrix)
clus <- hclust(distM, method = "average")
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
