library(cbioportalR)
library(dplyr)


## Link: https://www.karissawhiting.com/cbioportalR/articles/overview-of-workflow.html

set_cbioportal_db("public") #conectarnos a cbioportal
test_cbioportal_db() #testar que estamos conectados


all_studies <- available_studies() #descargarnos todos los estudios
all_studies


get_study_info("brca_tcga_pan_can_atlas_2018") %>%
  t() #Ver la información del estudio que nos estamos descargando

available_profiles(study_id = "brca_tcga_pan_can_atlas_2018") #ver tipo de archivos que hay

available_profiles(study_id = "brca_tcga_pan_can_atlas_2018") %>%
  pull(molecularProfileId) #ver archivos que hay en molecular profile


genetics <- get_genetics_by_study(study_id = "brca_tcga_pan_can_atlas_2018") # Get All Genomic Information By Study

attr_brca <- available_clinical_attributes("brca_tcga_pan_can_atlas_2018") #Get all available clinical attribute IDs for a study


p1 <- available_patients("brca_tcga_pan_can_atlas_2018") #Get All Patient IDs in a Study

clinical_brca <- get_clinical_by_patient(patient_id = p1$patientId,
                                         study_id = "brca_tcga_pan_can_atlas_2018")

head(clinical_brca)