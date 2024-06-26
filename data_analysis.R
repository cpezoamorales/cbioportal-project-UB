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
if(!require("RColorBrewer")){install.packages("RColorBrewer")}

if(!requireNamespace('BiocManager',quietly = T)) install.packages("BiocManager")

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


clinical_brca_wide$Fraction.Genome.Altered <- as.numeric (clinical_brca_wide$Fraction.Genome.Altered)


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

#Grupos de age
clinical_brca_wide <- clinical_brca_wide %>%
  mutate(age_group = case_when(
    AGE < 40 ~ "<40",
    AGE >= 40 & AGE <= 60 ~ "40-60",
    AGE > 60 ~ ">60"
  ))

# Ordenar los grupos
clinical_brca_wide$age_group <- factor(clinical_brca_wide$age_group,
                                       levels = c("<40", "40-60", ">60"))



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


## Analisis descriptivo en función de SUBTYPE

tbl_subtype <- tbl_summary(
  clinical_brca_wide %>%
    select(-c(patientId, FORM_COMPLETION_DATE, ICD_10, OTHER_PATIENT_ID )), 
  by = SUBTYPE, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()

tbl_subtype


## Analisis descriptivo en función de GROUP STAGE

tbl_groupstage <- tbl_summary(
  clinical_brca_wide %>%
    select(-c(patientId, FORM_COMPLETION_DATE, ICD_10, OTHER_PATIENT_ID )), 
  by = group_stage, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()

tbl_groupstage

## Analisis descriptivo en función de AGE GROUP
tbl_agegroup <- tbl_summary(
  clinical_brca_wide %>%
    select(-c(patientId, FORM_COMPLETION_DATE, ICD_10, OTHER_PATIENT_ID )), 
  by = age_group, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()

tbl_agegroup


## Analisis descriptivo en función de CancerType (pero solo los 2 grupos principales: Ductal y Lobular)

clinical_brca_wide_filter_cancertype <- clinical_brca_wide %>%
  filter(Cancer.Type.Detailed %in% c("Breast Invasive Ductal Carcinoma", "Breast Invasive Lobular Carcinoma")) #se crea subpoblacion con los 2 tipos de cancer más frecuentes para compararlos

tbl_cancertype <- tbl_summary(
  clinical_brca_wide_filter_cancertype %>%
    select(-c(patientId, FORM_COMPLETION_DATE, ICD_10, OTHER_PATIENT_ID)), 
  by = Cancer.Type.Detailed, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()

tbl_cancertype

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

#HR del modelo de cox
cox_pfs_subtype <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_subtype_tbl <- tbl_regression(cox_pfs_subtype, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_subtype_tbl


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


### PFS por subtipos  Zoom de los 3 primeros años (hasta 3 años = 36 meses)
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

#HR del modelo de cox
cox_pfs_stage <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_stage_tbl <- tbl_regression(cox_pfs_stage, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_stage_tbl




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

#HR del modelo de cox
cox_pfs_stage_group <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_stage_group_tbl <- tbl_regression(cox_pfs_stage_group, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_stage_group_tbl



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

#HR del modelo de cox
cox_pfs_lymphnodes <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_lymphnodes_tbl <- tbl_regression(cox_pfs_lymphnodes, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_lymphnodes_tbl

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

#HR del modelo de cox
cox_pfs_race <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ RACE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_race_tbl <- tbl_regression(cox_pfs_race, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_race_tbl

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


#HR del modelo de cox
cox_pfs_cancer_type <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ Cancer.Type.Detailed, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_cancer_type_tbl <- tbl_regression(cox_pfs_cancer_type, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_cancer_type_tbl


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


#HR del modelo de cox
cox_pfs_oncotree <- coxph(Surv(PFS_MONTHS, PFS_STATUS == "1") ~ Oncotree.Code, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_pfs_oncotree_tbl <- tbl_regression(cox_pfs_oncotree, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_pfs_oncotree_tbl





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

#HR del modelo de cox
cox_os_subtype <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_subtype_tbl <- tbl_regression(cox_os_subtype, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_subtype_tbl


### OS por subtipos  Zoom de los 10 primeros años (hasta 120 meses)
os_subtype_120m <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) +  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x



os_subtype_120m


### OS por subtipos  Zoom de los 3 primeros años (hasta 3 años = 36 meses)
os_subtype_36m <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) +  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 36))  # Ajustar los límites del eje x

os_subtype_36m




### OS por AJCC stage
os_stage <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 36))  # Ajustar los límites del eje x

os_stage

#HR del modelo de cox
cox_os_stage <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_stage_tbl <- tbl_regression(cox_os_stage, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_stage_tbl




### OS por AJCC stage en 120m
os_stage_120m <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x

os_stage_120m


### OS por Group stage stage
os_stage_group <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by GROUP STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 36))  # Ajustar los límites del eje x

os_stage_group

#HR del modelo de cox
cox_os_stage_group <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_stage_group_tbl <- tbl_regression(cox_os_stage_group, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_stage_group_tbl



### OS por Group stage stage en 120m
os_stage_120m_group <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by GROUP STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2) +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 120))  # Ajustar los límites del eje x

os_stage_120m_group


### OS por Lymph node
os_lympnodes <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by LYMPH NODES"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

os_lympnodes

#HR del modelo de cox
cox_os_lymphnodes <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_lymphnodes_tbl <- tbl_regression(cox_os_lymphnodes, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_lymphnodes_tbl

### OS por Race
os_race <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ RACE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by RACE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

os_race

#HR del modelo de cox
cox_os_race <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ RACE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_race_tbl <- tbl_regression(cox_os_race, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_race_tbl

### OS por Cancer Type
os_cancer_type <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ Cancer.Type.Detailed , data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by Cancer Type"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

os_cancer_type


#HR del modelo de cox
cox_os_cancer_type <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ Cancer.Type.Detailed, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_cancer_type_tbl <- tbl_regression(cox_os_cancer_type, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_cancer_type_tbl



### OS por Oncotree
os_oncotree <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ Oncotree.Code, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by CANCER TYPE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

os_oncotree


#HR del modelo de cox
cox_os_oncotree <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ Oncotree.Code, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_oncotree_tbl <- tbl_regression(cox_os_oncotree, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_oncotree_tbl


### OS por Age group
os_age <-survfit2(Surv(OS_MONTHS, OS_STATUS == "1") ~ age_group, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "OS by AGE GROUP"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

os_age

#HR del modelo de cox
cox_os_age <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ age_group, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_age_tbl <- tbl_regression(cox_os_age, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_age_tbl





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

# Crear data con solo filas que nos interesan
mutations_filtered <- mutations %>%
  select(hugoGeneSymbol, patientId, variantType)

# Calcular la prevalencia de mutaciones para cada gen y tipo de variante
gene_variant_prevalence <- mutations_filtered %>%
  group_by(hugoGeneSymbol, variantType) %>%
  summarise(prevalence = n(), .groups = 'drop') %>%
  arrange(desc(prevalence))

# Seleccionar las 50 mutaciones más prevalentes
top_mutations <- gene_variant_prevalence %>%
  slice_head(n = 50) %>%
  select(hugoGeneSymbol, variantType)


# Filtrar el data frame original para incluir solo las mutaciones seleccionadas
mutations_filtered_top <- mutations_filtered %>%
  semi_join(top_mutations, by = c("hugoGeneSymbol", "variantType"))

#Pivotar la tabla de larga a ancha
mutations_wide <- mutations_filtered_top %>%
  group_by(hugoGeneSymbol, patientId) %>%
  summarise(variantType = paste(unique(variantType), collapse = ";"), .groups = 'drop') %>%
  pivot_wider(names_from = patientId, values_from = variantType, values_fill = list(variantType = "WT"))

# Ver el resultado
print(mutations_wide)

# Crear una versión modificada de la tabla para el heatmap
mutations_wide_modified <- mutations_wide %>%
  pivot_longer(-hugoGeneSymbol, names_to = "patientId", values_to = "variantType") %>%
  mutate(variantType = factor(variantType, levels = c("WT", "SNP", "INS", "DEL", "CNV", "AMP", "HOMD"))) %>%
  pivot_wider(names_from = patientId, values_from = variantType)

# Convertir a matriz
mutation_matrix <- as.matrix(mutations_wide_modified %>% select(-hugoGeneSymbol))
mutation_matrix <- mutation_matrix[, order(colnames(mutation_matrix))]

# Asignar los nombres de las filas
rownames(mutation_matrix) <- mutations_wide_modified$hugoGeneSymbol

# Definir una paleta de colores personalizada
variant_colors <- c("WT" = "#FFF5EE", "SNP" = "#A2CD5A", "INS" = "#EE6363", "DEL" = "#6CA6CD")



# Crear el heatmap
Heatmap(mutation_matrix,
        name = "Mutations",
        col = variant_colors,
        show_row_names = TRUE,
        show_column_names = FALSE,
        column_title = "Patients",
        row_title = "Genes",
        cluster_rows = FALSE,    
        cluster_columns = FALSE, 
        show_heatmap_legend = TRUE, 
        show_column_dend = FALSE,  
        show_row_dend = FALSE,
)



###### Heatmap con la función específica de ONCOPRINT de cBioportal para resaltar también los porcentajes y ver el total de mutaciones y los tipos #########

# Definir la paleta de colores personalizada, excluyendo WT
col <- c("SNP" = "#A2CD5A", "INS" = "#EE6363", "DEL" = "#6CA6CD")

# Convertir valores WT a NA para que no los contabilize como mutaciones y se puedan calcular los % correctamente
mutation_matrix[mutation_matrix == "WT"] <- NA

# Crear una lista de funciones para las variantes
alter_fun <- list(
  background = alter_graphic("rect", fill = "#FFF5EE"),  # Fondo para WT
  SNP = alter_graphic("rect", fill = col["SNP"]),
  INS = alter_graphic("rect", fill = col["INS"]),
  DEL = alter_graphic("rect", fill = col["DEL"])
)

# Definir los parámetros de la leyenda del heatmap
heatmap_legend_param <- list(
  title = "Mutations",
  at = names(col),
  labels = c("SNP", "INS", "DEL")
)

# Crear el oncoPrint con genes en el eje Y y pacientes en el eje X
oncoPrint(mutation_matrix,
          alter_fun = alter_fun,
          col = col,
          top_annotation =
          show_row_names = TRUE,
          show_column_names = FALSE,
          column_title = "Patients",
          row_title = "Genes",
          heatmap_legend_param = heatmap_legend_param,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE
          )

##INTENTO DE PRUEBA CON DIVISON DE SUBTIPOS Y EDAD, FAIL

#Selecciono las variables de clinical info que me interesan y creo una base nueva
patient_info <- clinical_brca_wide %>%
  select(patientId, SUBTYPE, age_group)

patient_info$SUBTYPE<- as.factor(patient_info$SUBTYPE)
patient_info$age_group<- as.factor(patient_info$age_group)

# Definir colores para SUBTYPE y age_group
subtype_colors <- c("BRCA_Basal" = "#ff9634", "BRCA_Her2" = "#f99e61", "BRCA_LumA" = "#efa689", "BRCA_LumB" = "#e2aeb0", "BRCA_Normal" = "#cfb6d7")
age_group_colors <- c("<40" = "#FFF68F", ">60" = "#CDC673", "40-60" = "#8B864E")


# Función para generar barras de anotación basadas en SUBTYPE
anno_subtype <- function(patient_info) {
  HeatmapAnnotation(
    foo = anno_barplot(factor(patient_info$SUBTYPE), col = subtype_colors),
    width = unit(1, "cm")
  )
}

# Función para generar barras de anotación basadas en age_group
anno_age_group <- function(patient_info) {
  HeatmapAnnotation(
    foo = anno_barplot(factor(patient_info$age_group), col = age_group_colors),
    width = unit(1, "cm")
  )
}


# Definir la paleta de colores personalizada, excluyendo WT
col <- c("SNP" = "#A2CD5A", "INS" = "#EE6363", "DEL" = "#6CA6CD")

# Convertir valores WT a NA para que no los contabilize como mutaciones y se puedan calcular los % correctamente
mutation_matrix[mutation_matrix == "WT"] <- NA

# Crear la lista de funciones para las variantes (alter_fun)
alter_fun <- list(
  background = alter_graphic("rect", fill = "#FFF5EE"),  # Fondo para WT
  SNP = alter_graphic("rect", fill = col["SNP"]),
  INS = alter_graphic("rect", fill = col["INS"]),
  DEL = alter_graphic("rect", fill = col["DEL"])
)



# Crear el oncoPrint con genes en el eje Y y pacientes en el eje X
oncoPrint(
  mutation_matrix,
  top_annotation = rowAnnotation(
    SUBTYPE = anno_subtype(patient_info),
    age_group = anno_age_group(patient_info)
  ),
  alter_fun = alter_fun,
  col = col, 
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Patients",
  row_title = "Genes",
  heatmap_legend_param = list(
    title = "Mutations",
    at = c("SNP", "INS", "DEL"),
    labels = c("SNP", "Insertion", "Deletion")
  ),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE
)
