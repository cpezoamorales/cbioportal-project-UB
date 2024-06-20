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
if(!require("ggfortify")){install.packages("ggfortify")}
if(!require("gt")){install.packages("gt")}
if(!require("circlize")){install.packages("circlize")}
if(!require("patchwork")){install.packages("patchwork")}

if(!requireNamespace('BiocManager',quietly = T)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)



# Carga de datos crusdos almacenados previamente en archivo RData
load("data.RData")

data_oncotree <- read.csv2("Oncotree_Code.full.csv")
data_oncotree <- data_oncotree %>% rename(patientId = Patient.ID)

data_tmb <- read.csv2("TMB_(nonsynonymous).csv")
data_tmb <- data_tmb %>% rename(patientId = Patient.ID)

data_fraction_gen_al <- read.csv2("Mutation_Count_vs_Fraction_Genome_Altered.csv")
data_fraction_gen_al <- data_fraction_gen_al %>% rename(patientId = Patient.ID)

data_cancer_type <- read.csv2("Cancer_Type_Detailed.csv")
data_cancer_type <- data_cancer_type %>% rename(patientId = Patient.ID)

# Seleccionar solo las columnas necesarias de cada dataset para hacer el merge
data_oncotree <- data_oncotree %>% select(patientId, Oncotree.Code)
data_tmb <- data_tmb %>% select(patientId, TMB..nonsynonymous.)
data_fraction_gen_al <- data_fraction_gen_al %>% select(patientId, Fraction.Genome.Altered, Mutation.Count)
data_cancer_type <- data_cancer_type %>% select(patientId, Cancer.Type.Detailed)


# Transformar la tabla en Wide con pivot_wider, crea columnas viniendo de muchas filas
clinical_brca_wide <- clinical_brca %>%
  select(patientId, clinicalAttributeId, value)  %>%
  pivot_wider(names_from = clinicalAttributeId, values_from = value)

str(clinical_brca_wide)

length(unique(clinical_brca$patientId))


# Realizar los merges para añadir las variables Oncotree, TMB, Fraccion de gen alterada y cancer type que estaban en otras tablas a parte
clinical_brca_wide <- clinical_brca_wide %>%
  left_join(data_oncotree, by = "patientId") %>%
  left_join(data_tmb, by = "patientId") %>%
  left_join(data_fraction_gen_al, by = "patientId") %>%
  left_join(data_cancer_type, by = "patientId")




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



#### PFS #####

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
  add_risktable() + 
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
  add_risktable() + 
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
  add_risktable() + 
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

cna <- cna %>% 
  mutate(alteration_type = case_when(
    alteration == -2 ~ "HOMDEL",
    alteration == 2 ~ "AMP"
  ))



#Crear nova variable group Stage

clinical_brca_wide <- clinical_brca_wide %>% 
  mutate(group_stage = case_when(
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE I" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IA" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IB")  ~ "Stage I ",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE II" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIA" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIB")  ~ "Stage II ",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE III" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIA" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIB" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIC")  ~ "Stage III ",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IV") ~ "Stage IV ",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE X") ~ "Stage X",
))

print(as_tibble(clinical_brca_wide) %>%
        select(AJCC_PATHOLOGIC_TUMOR_STAGE, group_stage) %>%
        arrange(AJCC_PATHOLOGIC_TUMOR_STAGE), n = 1084)

### PFS por AJCC stage
pfs_stage <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 36))  # Ajustar los límites del eje x

pfs_stage



### PFS por AJCC stage en 120m
pfs_stage_120m <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x

pfs_stage_120m


### PFS por Group stage stage
pfs_stage_group <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by GROUP STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 36))  # Ajustar los límites del eje x

pfs_stage_group



### PFS por Group stage stage en 120m
pfs_stage_120m_group <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by GROUP STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x

pfs_stage_120m_group


### PFS por Lymph node
pfs_lympnodes <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by LYMPH NODES"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

pfs_lympnodes

### PFS por Race
pfs_race <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ RACE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by RACE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

pfs_race

### PFS por Cancer Type
pfs_cancer_type <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ Cancer.Type.Detailed , data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by Cancer Type"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

pfs_cancer_type


clinical_brca_wide$Oncotree.Code <- as.factor(clinical_brca_wide$Oncotree.Code)

### PFS por Oncotree
pfs_oncotree <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ Oncotree.Code, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by CANCER TYPE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

pfs_oncotree


#Grupos de age
clinical_brca_wide <- clinical_brca_wide %>%
  mutate(age_group = case_when(
    AGE < 40 ~ "<40",
    AGE >= 40 & AGE <= 60 ~ "40-60",
    AGE > 60 ~ ">60"
  ))


### PFS por Age group
pfs_age <-survfit2(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ age_group, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by AGE GROUP"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

pfs_age

#HR del modelo de cox
cox_pfs_age <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ age_group, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_age_tbl <- tbl_regression(cox_pfs_age, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_age_tbl



combined_plot <- pfs_age + as_gt(cox_pfs_age_tbl) + plot_layout(ncol = 1)


#### OS #####

# Cambiar los niveles de la variable OS_STATUS
clinical_brca_wide$OS_STATUS <- factor(clinical_brca_wide$OS_STATUS, 
                                        levels = c("0:LIVING", "1:DECEASED"), 
                                        labels = c("0", "1"))

# Verificar los nuevos niveles
levels(clinical_brca_wide$OS_STATUS)


#ALTERNATIVA paquete  "ggsurvfit" --> mejor este se puede editar con ggplot2
os <- survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ 1, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1, color = "#607B8B") +  # Cambiar el color de la línea
  add_confidence_interval(fill = "#B0E2FF") +  # Cambiar el color del sombreado que representa el CI
  add_risktable() +  # Cambiar el estilo de la tabla de riesgo
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Overall survival",
    x= " Time (months)",
    title = "OS"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) 

os


### OS por subtipos
os_subtype <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Overall survival",
    x= " Time (months)",
    title = "OS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2)  # Añadir marcas de censura

os_subtype



####### VARIABLES NUMERICAS ###########

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

# MUTATIONS analysis
mutations$mutationType <- as.factor(mutations$mutationType)
levels(mutations$mutationType)
mutations$variantType <- as.factor(mutations$variantType)
levels(mutations$variantType)

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

table_mutations <- summary_data %>%
  gt() %>%
  cols_hide(columns = "total_samples")

table_mutations

# STRUCTURAL VAARIANTS analysis
summary_data_structural <- structural_variant %>%
  group_by(site1HugoSymbol) %>%
  summarise(
    total_mutations = n(),  # Número total de SV en este gen
    unique_samples = n_distinct(sampleId),  # Número de muestras únicas con al menos una SV en este gen
    total_samples = n_distinct(structural_variant$sampleId)  # Número total de muestras en todo el dataset
  ) %>%
  mutate(
    percent_samples_with_sv = (unique_samples / total_samples) * 100  # Porcentaje de muestras con al menos una SV  en este gen
  ) %>%
  arrange(desc(percent_samples_with_sv)) %>%
  mutate(percent_samples_with_sv = round(percent_samples_with_sv, 2))


table_sv <- summary_data_structural %>%
  gt() %>%
  cols_hide(columns = "total_samples")

table_sv

#COPY NUMBER ALTERATIONS

cna$alteration <- as.factor(cna$alteration)

levels(cna$alteration)

cna <- cna %>% 
  mutate(alteration_type = case_when(
    alteration == -2 ~ "HOMDEL",
    alteration == 2 ~ "AMP"
  ))


summary_data_cna <- cna %>% #HAY UN PROBLEMA QUE LO HACE PERO TARDA MUCHO EN CALCULARLO porque hay >446.000 lineas
  group_by(hugoGeneSymbol, alteration_type) %>%
  summarise(
    total_mutations = n(),  # Número total de CNA en este gen
    unique_samples = n_distinct(sampleId),  # Número de muestras únicas con al menos una CNA en este gen
    total_samples = n_distinct(cna$sampleId)
    # Número total de muestras en todo el dataset
  ) %>%
  mutate(
    percent_samples_with_cna = (unique_samples / total_samples) * 100  # Porcentaje de muestras con al menos una CNA  en este gen
  ) %>%
  arrange(desc(percent_samples_with_cna)) %>%
  mutate(percent_samples_with_cna = round(percent_samples_with_cna, 2))


table_cna <- summary_data_cna %>%
  gt() %>%
  cols_hide(columns =  c("total_samples", "total_mutations"))

table_cna


############### HEATMAP ####################
#Provar de hacer un heatmap con mutaciones (mutationType i variantType --> mirar les variables que son una N)

structural_variant$site2EffectOnFrame <- as.factor(structural_variant$site2EffectOnFrame)
