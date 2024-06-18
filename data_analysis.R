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


### Análisis exploratorios: Correlación de variables numéricas
clinical_brca_numeric <- clinical_brca_wide %>% 
  select_if(is.numeric)

# lower case:
names(clinical_brca_numeric) <- c("age", "buffa_hypoxia_score", "days_last_followup",
                                  "days_birth", "days_pathologic_diag", "dfs_month",
                                  "dss_month", "os_month", "pfs_month", 
                                  "ragnum_hypoxia_score", "winter_hypoxia_score")

pairs(clinical_brca_numeric)
plot(clinical_brca_numeric$days_pathologic_diag)

clinical_brca_numeric_2 <- clinical_brca_numeric
clinical_brca_numeric_2$days_pathologic_diag <- NULL
pairs(clinical_brca_numeric_2)

clinical_corr <- cor(clinical_brca_numeric_2)
corrplot(clinical_corr, type = "upper")
# problema al hacer correlación, quizás por distinto rango de los valores
# Normalizar los datos?

par(mfrow=c(1,3))
boxplot(clinical_brca_numeric$age)
hist(clinical_brca_numeric$age)
plot(clinical_brca_numeric$age, pch=20)

boxplot(clinical_brca_numeric$winter_hypoxia_score)
hist(clinical_brca_numeric$winter_hypoxia_score)
hist(clinical_brca_numeric$ragnum_hypoxia_score)
hist(clinical_brca_numeric$buffa_hypoxia_score)
hist(clinical_brca_numeric$dfs_month)

###--------------------------------------------------------------------------###
### Analisis de expresión de genes con WGCNA
###--------------------------------------------------------------------------###
# Instalar Bioconductor
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.19")

# Instalar y cargar el paquete WGCNA
BiocManager::install("WGCNA")
library("WGCNA")

# Lectura de datos
rna_data <- read_csv2("data_mrna_seq_v2_rsem.csv")
col_select <- rna_data[, -1]

# Transformación de tabla para el análisis
rna_longer <- t(col_select)
colnames(rna_longer) <- rna_longer[1,]
rna_longer <- rna_longer[-1, ]

# look for missing values
gsg = goodSamplesGenes(rna_longer, verbose = 3);
gsg$allOk

sampleTree <- hclust(dist(rna_longer), method = "average")
plot(sampleTree)

clust = cutreeStatic(sampleTree, cutHeight = 31000, minSize = 10)
