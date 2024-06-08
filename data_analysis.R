# Load data, this format is faster thar read from the cbioportal API
load(data.RData)

# Transform the table clinical_brca:
data_1 <- clinical_brca %>% 
  select(patientId, clinicalAttributeId, value) %>%
  pivot_wider(names_from = clinicalAttributeId, values_from = value)

str(data_1)
# Transform numerical values from chr to int