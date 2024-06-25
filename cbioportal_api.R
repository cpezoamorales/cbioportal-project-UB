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

# Guardar los datos de interés para cargar rápidamente en una futura sesión de R:
save(genetics, clinical_brca, attr_brca, file = "data.RData")
