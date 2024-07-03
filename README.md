# cbioportal-project-UB
Repositorio para el trabajo final del grupo 1 del curso de Ciencia de Datos del Instituto de Formación Continua-IL3 Universidad de Barcelona.

El objetivo de este trabajo fue aplicar técnicas de análisis de datos en un set de datos biomédicos obtenidos desde el portal cBioPortal, un portal que almacena los datos de varios estudios de cancer a través de la colaboración de distintos grupos de investigación a nivel mundial. 

El primer paso consistió en la búsqueda de un dataset en cBioportal con las características que nos interesan: multiples variables y número de pacientes superior a 500. El estudio de interés sobre cancer de mama de interés se encuentra en el siguiente link: https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018.

Este repositorio contiene los códigos en R generados en este trabajo los cuales son los siguientes:
1) cbioportal_api.R = código que permite conectar una sesión de R con la API web del portal cbioportal permitiendo la descarga de los datos disponible. Además incluye un preprocesamiento de los datos y almacenamiento en formato RData de todas las tablas creadas para su posterior análisis. 
2) cbioportal_data.RData = archivo que contiene los datos de interés para desarrollar los análisis.
3) data_analysis.R = código que contiene el análisis descirptivo de los datos, así como tambien análisis de curvas se sobrevivencia, análisis con métodos de aprendizaje automático no supervisado y oncoplots.
4) app_survival.R = código de la aplicación de Shiny para la visualización interactiva de las curvas de superviviencia.

Con estos códigos se aborda metodológicamente los objetivos específicos del proyecto:
1) Conectar cBioPortal con una sesión de  R mediante una interfaz-API para obtención de datos.
2) Realizar estadística descriptiva de los datos clínicos generando tablas y gráficos para una fácil visualización.
3) Generar curvas de supervivencia Kaplan-Meier en la población general y en subgrupos específicos.
4) Buscar patrones en los datos mediante técnicas de correlación y aprendizaje automático no supervisado.
5) Visualizar patrones de alteraciones genéticas prevalentes en la población de estudio mediante un Heatmap.
6) Desarrollar una interfaz interactiva utilizando Shiny app para la visualización de una parte de los datos. 
