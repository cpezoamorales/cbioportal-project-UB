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

# Carga de datos crusdos almacenados previamente en archivo RData para que no tarde tanto en cargar los datos
load("data2.RData")


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
  italicize_levels() %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#EE6363'>Análisis exploratorio de variables</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )


tbl

# Filtrando las variables numericas seleccionadas
selected_vars_num <- c("AGE", "BUFFA_HYPOXIA_SCORE", "RAGNUM_HYPOXIA_SCORE", "WINTER_HYPOXIA_SCORE",
  "TMB..nonsynonymous.", "Fraction.Genome.Altered", "Mutation.Count")

tbl_selection_num <- tbl_summary(
  clinical_brca_wide %>%
    select(all_of(selected_vars_num)), 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_header(label = "**Numeric characteristics**") %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#7AC5CD'>Análisis exploratorio de variables numéricas</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )

tbl_selection_num


# Filtrando las variables categóricas seleccionadas
selected_vars_cat <- c("SEX", "age_group", "ETHNICITY", "RACE", "SUBTYPE", "Cancer.Type.Detailed",
  "AJCC_PATHOLOGIC_TUMOR_STAGE", "group_stage", "PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT",
  "PATH_M_STAGE", "PRIOR_DX", "HISTORY_NEOADJUVANT_TRTYN",
  "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT", "DFS_STATUS", "PFS_STATUS",
  "DSS_STATUS", "OS_STATUS")

tbl_selection_cat <- tbl_summary(
  clinical_brca_wide %>%
    select(all_of(selected_vars_cat)), 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels()  %>%
  modify_header(label = "**Categorical characteristics**") %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#FF7F00'>Análisis exploratorio de variables categóricas</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )

tbl_selection_cat


## Analisis descriptivo en función de SUBTYPE

tbl_subtype <- tbl_summary(
  clinical_brca_wide %>%
    select(all_of(c(selected_vars_num, selected_vars_cat))), 
  by = SUBTYPE, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()  %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#A2CD5A'>Descriptiva por Subtipo</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )

tbl_subtype


## Analisis descriptivo en función de GROUP STAGE

tbl_groupstage <- tbl_summary(
  clinical_brca_wide %>%
    select(all_of(c(selected_vars_num, selected_vars_cat))), 
  by = group_stage, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()  %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#A2CD5A'>Descriptiva por Grupo de Estadío</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )

tbl_groupstage

## Analisis descriptivo en función de AGE GROUP
tbl_agegroup <- tbl_summary(
  clinical_brca_wide %>%
    select(all_of(c(selected_vars_num, selected_vars_cat))), 
  by = age_group, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()  %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#A2CD5A'>Descriptiva por Grupo de Edad</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )

tbl_agegroup


## Analisis descriptivo en función de CancerType (pero solo los 2 grupos principales: Ductal y Lobular)

clinical_brca_wide_filter_cancertype <- clinical_brca_wide %>%
  filter(Cancer.Type.Detailed %in% c("Breast Invasive Ductal Carcinoma", "Breast Invasive Lobular Carcinoma")) #se crea subpoblacion con los 2 tipos de cancer más frecuentes para compararlos

tbl_cancertype <- tbl_summary(
  clinical_brca_wide_filter_cancertype %>%
    select(all_of(c(selected_vars_num, selected_vars_cat))), 
  by = Cancer.Type.Detailed, 
  type = list(where(is.factor) ~ "categorical"),
  missing_text = "(Missing)",
  sort = list(everything() ~ "alphanumeric"),
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{median}, Range ({min} - {max})")
) %>%
  add_n() %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_p()  %>%
  as_gt() %>%
  tab_header(
    title = md("<span style='font-size:24px; color:#A2CD5A'>Descriptiva por Tipo Histológico</span>")
  ) %>%
  tab_options(
    heading.title.font.size = px(24), 
    heading.title.font.weight = "bold"
  )


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
  coord_cartesian(xlim = c(0, 200))  # Ajustar los límites del eje x

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
  add_risktable() + ##Quitarlo cuando no lo queramos
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Progression free survival",
    x= " Time (months)",
    title = "PFS by GROUP STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)  +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 200))  # Ajustar los límites del eje x


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
  add_risktable() + 
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
    y = "Overall survival",
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
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Overall survival",
    x= " Time (months)",
    title = "OS by AGE GROUP"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

os_age

#HR del modelo de cox
cox_os_age <- coxph(Surv(OS_MONTHS, OS_STATUS == "1") ~ age_group, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_os_age_tbl <- tbl_regression(cox_os_age, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_os_age_tbl

#### DSS #####

# Cambiar los niveles de la variable DSS_STATUS
clinical_brca_wide$DSS_STATUS <- factor(clinical_brca_wide$DSS_STATUS, 
                                        levels = c("0:ALIVE OR DEAD TUMOR FREE", "1:DEAD WITH TUMOR"), 
                                        labels = c("0", "1"))

# Verificar los nuevos niveles
levels(clinical_brca_wide$DSS_STATUS)

#ALTERNATIVA paquete  "ggsurvfit" --> mejor este se puede editar con ggplot2
dss <- survfit2(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ 1, data = clinical_brca_wide) %>%
  ggsurvfit(linewidth = 1, color = "#607B8B") +  # Cambiar el color de la línea
  add_confidence_interval(fill = "#B0E2FF") +  # Cambiar el color del sombreado que representa el CI
  add_risktable() +  # Cambiar el estilo de la tabla de riesgo
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Disease Specific survival",
    x= " Time (months)",
    title = "DSS"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2) 

dss




### DSS por subtipos
dss_subtype <-survfit2(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Disease Specific survival",
    x= " Time (months)",
    title = "DSS by subtype"
  ) +
  add_censor_mark(color = "#607B8B", shape = 124, size = 2)  # Añadir marcas de censura

dss_subtype

#HR del modelo de cox
cox_dss_subtype <- coxph(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ SUBTYPE, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_dss_subtype_tbl <- tbl_regression(cox_dss_subtype, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_dss_subtype_tbl


### DSS por Group stage stage
dss_stage_group <-survfit2(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Disease Specific survival",
    x= " Time (months)",
    title = "DSS by GROUP STAGE"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)  +  # Añadir marcas de censura  # Añadir marcas de censura
  coord_cartesian(xlim = c(0, 200))  # Ajustar los límites del eje x


dss_stage_group

#HR del modelo de cox
cox_dss_stage_group <- coxph(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ group_stage, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_dss_stage_group_tbl <- tbl_regression(cox_dss_stage_group, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_dss_stage_group_tbl

### DSS por Age group
dss_age <-survfit2(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ age_group, data = clinical_brca_wide) %>% 
  ggsurvfit(linewidth = 1)  + 
  add_risktable() + 
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() +
  labs(
    y = "Disease Specific survival",
    x= " Time (months)",
    title = "DSS by AGE GROUP"
  ) +
  add_censor_mark(color = "gray50", shape = 124, size = 2)

dss_age

#HR del modelo de cox
cox_dss_age <- coxph(Surv(DSS_MONTHS, DSS_STATUS == "1") ~ age_group, data = clinical_brca_wide) # Ajustar el modelo de Cox

cox_dss_age_tbl <- tbl_regression(cox_dss_age, exponentiate = TRUE) # Crear una tabla de resumen del modelo de Cox usando tbl_regression

cox_dss_age_tbl





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
par(mfrow=c(2,5)) #para determinar estructura gráficos 
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
  mutate(percent_samples_with_mutation = round(percent_samples_with_mutation, 2))%>%
  slice_head(n = 100)  # Seleccionar los top 100

table_mutations <- summary_data %>%
  gt() %>%
  cols_hide(columns = "total_samples") %>%
  tab_header(
    title = md("<span style='font-size:20px; color:#EE6363'>Resumen de Mutaciones por Gen - Top 100 </span>")
  ) %>%
  cols_label(
    hugoGeneSymbol = "Gen",
    total_mutations = "Total de Mutaciones",
    unique_samples = "Total de Mutaciones por muestra única",
    percent_samples_with_mutation = "Porcentaje de Muestras con Mutación en este gen"
  )%>%
  fmt_number(
    columns = vars(percent_samples_with_mutation),
    decimals = 2,
    suffixing = TRUE
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = vars(hugoGeneSymbol))
  ) %>%
  fmt(
    columns = vars(percent_samples_with_mutation),
    fns = function(x) gsub("\\.", ",", paste0(x, "%"))
  )


table_mutations

# STRUCTURAL VARIANTS analysis
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
  mutate(percent_samples_with_sv = round(percent_samples_with_sv, 2)) %>%
  slice_head(n = 100)  # Seleccionar los top 100


table_sv <- summary_data_structural %>%
  gt() %>%
  cols_hide(columns = "total_samples")  %>%
  tab_header(
    title = md("<span style='font-size:20px; color:#EE6363'>Resumen de Variantes Estructurales por Gen - Top 100</span>")
  ) %>%
  cols_label(
    site1HugoSymbol	 = "Gen",
    total_mutations = "Total de Variantes estructurales",
    unique_samples = "Total de Variantes estructurales por muestra única",
    percent_samples_with_sv = "Porcentaje de Muestras con Variantes estructurales en este gen"
  )%>%
  fmt_number(
    columns = vars(percent_samples_with_sv),
    decimals = 2,
    suffixing = TRUE
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = vars(
      site1HugoSymbol))
  ) %>%
  fmt(
    columns = vars(percent_samples_with_sv),
    fns = function(x) gsub("\\.", ",", paste0(x, "%"))
  )


table_sv

#COPY NUMBER ALTERATIONS

cna$alteration <- as.factor(cna$alteration)

levels(cna$alteration)

cna <- cna %>% 
  mutate(alteration_type = case_when(
    alteration == -2 ~ "HOMDEL",
    alteration == 2 ~ "AMP",
    TRUE ~ NA_character_  # Agregar TRUE ~ NA_character_ para manejar otros valores posibles de alteration
  ))


summary_data_cna <- cna %>% 
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
  mutate(percent_samples_with_cna = round(percent_samples_with_cna, 2)) %>%
  ungroup() %>%  # Desagrupar antes de seleccionar los top 100
  slice_head(n = 100)  # Seleccionar los top 100


table_cna <- summary_data_cna %>%
  gt() %>%
  cols_hide(columns =  c("total_samples", "total_mutations")) %>%
  tab_header(
    title = md("<span style='font-size:20px; color:#EE6363'>Resumen de Alteraciones de número de copia (CNA) - Top 100</span>")
  ) %>%
  cols_label(
    hugoGeneSymbol	 = "Gen",
    total_mutations = "Total de Alteraciones de número de copia (CNA)",
    unique_samples = "Total de Alteraciones de número de copia (CNA) por muestra única",
    percent_samples_with_cna = "Porcentaje de Muestras con Alteraciones de número de copia (CNA) en este gen"
  )%>%
  fmt_number(
    columns = vars(percent_samples_with_cna),
    decimals = 2,
    suffixing = TRUE
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = vars(
      hugoGeneSymbol))
  ) %>%
  fmt(
    columns = vars(percent_samples_with_cna),
    fns = function(x) gsub("\\.", ",", paste0(x, "%"))
  )


table_cna


############### HEATMAP ####################

#Prueba Heatmap tipo 1 
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
          show_row_names = TRUE,
          show_column_names = FALSE,
          column_title = "Patients",
          row_title = "Genes",
          heatmap_legend_param = heatmap_legend_param,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE
)



#HEATMAP POR MUTATION TYPE

#Prueba Heatmap tipo 1
# Crear data con solo filas que nos interesan
mutations_filtered_2 <- mutations %>%
  select(hugoGeneSymbol, patientId, mutationType)

# Calcular la prevalencia de mutaciones para cada gen y tipo de variante
gene_variant_prevalence_2 <- mutations_filtered_2 %>%
  group_by(hugoGeneSymbol, mutationType) %>%
  summarise(prevalence = n(), .groups = 'drop') %>%
  arrange(desc(prevalence))

# Seleccionar las 50 mutaciones más prevalentes
top_mutations_2 <- gene_variant_prevalence_2 %>%
  slice_head(n = 50) %>%
  select(hugoGeneSymbol, mutationType)


# Filtrar el data frame original para incluir solo las mutaciones seleccionadas
mutations_filtered_2_top <- mutations_filtered_2 %>%
  semi_join(top_mutations_2, by = c("hugoGeneSymbol", "mutationType"))

#Pivotar la tabla de larga a ancha
mutations_wide_2 <- mutations_filtered_2_top %>%
  group_by(hugoGeneSymbol, patientId) %>%
  summarise(mutationType = paste(unique(mutationType), collapse = ";"), .groups = 'drop') %>%
  pivot_wider(names_from = patientId, values_from = mutationType, values_fill = list(mutationType = "WT"))

# Ver el resultado
print(mutations_wide_2)

# Crear una versión modificada de la tabla para el heatmap
mutations_wide_2_modified <- mutations_wide_2 %>%
  pivot_longer(-hugoGeneSymbol, names_to = "patientId", values_to = "mutationType") %>%
  mutate(mutationType = factor(mutationType, levels = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site","Translation_Start_Site"
  ))) %>%
  pivot_wider(names_from = patientId, values_from = mutationType)

# Convertir a matriz
mutation_matrix_2 <- as.matrix(mutations_wide_2_modified %>% select(-hugoGeneSymbol))
mutation_matrix_2 <- mutation_matrix_2[, order(colnames(mutation_matrix_2))]

# Asignar los nombres de las filas
rownames(mutation_matrix_2) <- mutations_wide_2_modified$hugoGeneSymbol

# Definir una paleta de colores personalizada
variant_colors_2 <- c("Frame_Shift_Del" = "#FFA07A", "Frame_Shift_Ins" = "#79a2ab", "In_Frame_Del" = "#E49183", "In_Frame_Ins" = "#E8A49A", "Missense_Mutation" = "#b3577f", "Nonsense_Mutation" = "#CDB79E", "Splice_Region" = "#D35C79", "Splice_Site" = "#CE9F51", "Translation_Start_Site" = "#46307E" )



# Crear el heatmap
Heatmap(mutation_matrix_2,
        name = "Mutations",
        col = variant_colors_2,
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
col_2 <- c("Frame_Shift_Del" = "#FFA07A", "Frame_Shift_Ins" = "#79a2ab", "In_Frame_Del" = "#E49183", "In_Frame_Ins" = "#E8A49A", "Missense_Mutation" = "#b3577f", "Nonsense_Mutation" = "#CDB79E", "Splice_Region" = "#D35C79", "Splice_Site" = "#CE9F51", "Translation_Start_Site" = "#46307E" )

# Convertir valores WT a NA para que no los contabilice como mutaciones y se puedan calcular los % correctamente
mutation_matrix_2[mutation_matrix_2 == "WT"] <- NA

# Convertir valores NA a cadenas vacías
mutation_matrix_2[is.na(mutation_matrix_2)] <- ""

# Crear una lista de funciones para las variantes
alter_fun_2 <- list(
  background = alter_graphic("rect", fill = "#FFF5EE"),  # Fondo para WT
  Frame_Shift_Del = alter_graphic("rect", fill = col_2["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect", fill = col_2["Frame_Shift_Ins"]),
  In_Frame_Del = alter_graphic("rect", fill = col_2["In_Frame_Del"]),
  In_Frame_Ins = alter_graphic("rect", fill = col_2["In_Frame_Ins"]),
  Missense_Mutation = alter_graphic("rect", fill = col_2["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect", fill = col_2["Nonsense_Mutation"]),
  Splice_Region = alter_graphic("rect", fill = col_2["Splice_Region"]),
  Splice_Site = alter_graphic("rect", fill = col_2["Splice_Site"]),
  Translation_Start_Site = alter_graphic("rect", fill = col_2["Translation_Start_Site"])
)

# Definir los parámetros de la leyenda del heatmap
heatmap_legend_param <- list(
  title = "Mutations",
  at = names(col_2),
  labels = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Region", "Splice_Site", "Translation_Start_Site")
)

# Crear el oncoPrint con genes en el eje Y y pacientes en el eje X
oncoPrint(
  mutation_matrix_2,
  alter_fun = alter_fun_2,
  col = col_2,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Patients",
  row_title = "Genes",
  heatmap_legend_param = heatmap_legend_param,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE
)

### FAIL ### INTENTO DE PRUEBA CON DIVISON DE SUBTIPOS Y EDAD

#Selecciono las variables de clinical info que me interesan y creo una base nueva
patient_info <- clinical_brca_wide %>%
  select(patientId, SUBTYPE, age_group)

patient_info$SUBTYPE<- as.factor(patient_info$SUBTYPE)
patient_info$age_group<- as.factor(patient_info$age_group)

# Definir colores para SUBTYPE y age_group
subtype_colors <- c("BRCA_Basal" = "#ff9634", "BRCA_Her2" = "#f99e61", "BRCA_LumA" = "#efa689", "BRCA_LumB" = "#e2aeb0", "BRCA_Normal" = "#cfb6d7")
age_group_colors <- c("<40 years old" = "#FFF68F", "40-60 years old" = "#CDC673", ">60 years old" = "#8B864E")


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
