###--------------------------------------------------------------------------###
### Código para conexión de la sesión de R con la web API de cBioPortal
###--------------------------------------------------------------------------###

### Set directorio para trabajo local o usar una proyecto de RStudio
# Cambiar la frase "/dirección/de/trabajo/" con la dirección local para ejecución del código:
setwd("/dirección/de/trabajo/")
getwd()

# Instalación y carga de las librerias necesarias para este código:
if(!require("cbioportalR")){install.packages("cbioportalR")}
if(!require("tidyverse")){install.packages("tidyverse")}

# Conectar con el portal
set_cbioportal_db("public")
# Evaluar conexión:
test_cbioportal_db()

# Llamar y visualizar todos los estudios diponibles para descargar:
all_studies <- available_studies()
all_studies

# Obtener información del data set de interés
get_study_info("brca_tcga_pan_can_atlas_2018") %>%
  t()

# Perfiles de alteraciones de genes, de expresión de genes y proteínas contenidas en el estudio data set:
available_profiles(study_id = "brca_tcga_pan_can_atlas_2018")

available_profiles(study_id = "brca_tcga_pan_can_atlas_2018") %>%
  pull(molecularProfileId)

# Crear una objeto tipo lista con los datos genéticos del data set, descargámnolos desde la web:
genetics <- get_genetics_by_study(study_id = "brca_tcga_pan_can_atlas_2018")
attr_brca <- available_clinical_attributes("brca_tcga_pan_can_atlas_2018")

# Obtención del Id de los pacientes del dataset:
p1 <- available_patients("brca_tcga_pan_can_atlas_2018")

# Obtención de los datos clínicos de los pacientes:
clinical_brca <- get_clinical_by_patient(patient_id = p1$patientId,
                                         study_id = "brca_tcga_pan_can_atlas_2018")

#Obtener archivos con otras variables guardados en tablas
### Laia subió los datos en formato csv 
data_oncotree <- read.csv2("Oncotree_Code.full.csv")
data_oncotree <- data_oncotree %>% rename(patientId = Patient.ID)

data_tmb <- read.csv2("TMB_(nonsynonymous).csv")
data_tmb <- data_tmb %>% rename(patientId = Patient.ID)

data_fraction_gen_al <- read.csv2("Mutation_Count_vs_Fraction_Genome_Altered.csv")
data_fraction_gen_al <- data_fraction_gen_al %>% rename(patientId = Patient.ID)

data_cancer_type <- read.csv2("Cancer_Type_Detailed.csv")
data_cancer_type <- data_cancer_type %>% rename(patientId = Patient.ID)

###--------------------------------------------------------------------------###
### Curacion de datos
###--------------------------------------------------------------------------###

# Seleccionar solo las columnas necesarias de cada dataset para hacer el merge
data_oncotree <- data_oncotree %>% 
  select(patientId, Oncotree.Code)

data_tmb <- data_tmb %>% 
  select(patientId, TMB..nonsynonymous.)

data_fraction_gen_al <- data_fraction_gen_al %>% 
  select(patientId, Fraction.Genome.Altered, Mutation.Count)

data_cancer_type <- data_cancer_type %>% 
  select(patientId, Cancer.Type.Detailed)

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

# Crear nova variable group Stage
clinical_brca_wide <- clinical_brca_wide %>% 
  mutate(group_stage = case_when(
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE I" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IA" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IB")  ~ "Stage I",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE II" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIA" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIB")  ~ "Stage II",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE III" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIA" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIB" | AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IIIC")  ~ "Stage III",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IV") ~ "Stage IV",
    (AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE X") ~ "Stage X",
  ))

print(as_tibble(clinical_brca_wide) %>%
        select(AJCC_PATHOLOGIC_TUMOR_STAGE, group_stage) %>%
        arrange(AJCC_PATHOLOGIC_TUMOR_STAGE), n = 1084)

# Grupos de age
clinical_brca_wide <- clinical_brca_wide %>%
  mutate(age_group = case_when(
    AGE < 40 ~ "<40 years old",
    AGE >= 40 & AGE <= 60 ~ "40-60 years old",
    AGE > 60 ~ ">60 years old"
  ))

# Ordenar los grupos
clinical_brca_wide$age_group <- factor(clinical_brca_wide$age_group,
                                       levels = c("<40 years old", "40-60 years old", ">60 years old"))

# Guardar los datos de interés para cargar rápidamente en una futura sesión de R:
save(genetics, 
     clinical_brca, 
     attr_brca, 
     data_oncotree, 
     data_tmb, 
     data_fraction_gen_al, 
     data_cancer_type, 
     clinical_brca_wide, 
     file = "cbioportal_data.RData")
