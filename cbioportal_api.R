library(cbioportalR)
library(dplyr)


## Link: https://www.karissawhiting.com/cbioportalR/articles/overview-of-workflow.html

set_cbioportal_db("public")
test_cbioportal_db()


all_studies <- available_studies()
all_studies


get_study_info("brca_tcga_pan_can_atlas_2018") %>%
  t()

available_profiles(study_id = "brca_tcga_pan_can_atlas_2018")

available_profiles(study_id = "brca_tcga_pan_can_atlas_2018") %>%
  pull(molecularProfileId)


genetics <- get_genetics_by_study(study_id = "brca_tcga_pan_can_atlas_2018")

attr_brca <- available_clinical_attributes("brca_tcga_pan_can_atlas_2018")


p1 <- available_patients("brca_tcga_pan_can_atlas_2018")

clinical_brca <- get_clinical_by_patient(patient_id = p1$patientId,
                                         study_id = "brca_tcga_pan_can_atlas_2018")

