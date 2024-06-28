###--------------------------------------------------------------------------###
### Código para conexión de la sesión de R con la web API de cBioPortal
###--------------------------------------------------------------------------###

# Instalación y carga de las librerias necesarias para este código:
if(!require("cbioportalR")){install.packages("cbioportalR")}
if(!require("dplyr")){install.packages("dplyr")}

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

# Crear una objeto tipo lista con los datos genéticos del data set, descargámdolos desde la web:
genetics <- get_genetics_by_study(study_id = "brca_tcga_pan_can_atlas_2018")

attr_brca <- available_clinical_attributes("brca_tcga_pan_can_atlas_2018")

# Obtención del Id de los pacientes del dataset:
p1 <- available_patients("brca_tcga_pan_can_atlas_2018")

# Obtención de los datos clínicos de los pacientes:
clinical_brca <- get_clinical_by_patient(patient_id = p1$patientId,
                                         study_id = "brca_tcga_pan_can_atlas_2018")

#Obtener archivos con otras variables guardados en tablas

data_oncotree <- read.csv2("Oncotree_Code.full.csv")
data_oncotree <- data_oncotree %>% rename(patientId = Patient.ID)

data_tmb <- read.csv2("TMB_(nonsynonymous).csv")
data_tmb <- data_tmb %>% rename(patientId = Patient.ID)

data_fraction_gen_al <- read.csv2("Mutation_Count_vs_Fraction_Genome_Altered.csv")
data_fraction_gen_al <- data_fraction_gen_al %>% rename(patientId = Patient.ID)

data_cancer_type <- read.csv2("Cancer_Type_Detailed.csv")
data_cancer_type <- data_cancer_type %>% rename(patientId = Patient.ID)


# Guardar los datos de interés para cargar rápidamente en una futura sesión de R:
save(genetics, clinical_brca, attr_brca, data_oncotree, data_tmb, data_fraction_gen_al, data_cancer_type, file = "data2.RData")
