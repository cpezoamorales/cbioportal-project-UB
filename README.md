# cbioportal-project-UB
Repositorio para el trabajo final del grupo 1 del curso de Ciencia de Datos del Instituto de Formación Continua-IL3 Universidad de Barcelona.

El objetivo de este trabajo fue aplicar técnicas de análisis de datos en datos biomédicos obtenidos desde el portal cBioPortal, un portal que alamcena los datos de distintos estudios de cancer almacenados a través de la colaboración de distintos grupo de investigación a nivel mundial. 

Este repositorio contiene los scripts generados en este trabajo, los cuales incliuyen los siguientes pasos:
1) Búsqueda de un dataset en cBioportal con las características que nos interesan (multiples variables y número de pacientes superior a 500). El estudio de interés sobre cancer de mama de interés se encuentra en el siguiente link: https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018
2) Conectar cBioPortal con la sesión de RStudio mediante una interfaz-API para obtener el dataset escogido.
3) Exploración de los datos con R:
   3.1) Generación de tablas con estadísticas descriptivas para la mayoría de variables, con el cálculo de medidas de tendencia central para variables como: edad al diagnóstico, alteraciones genómicas, subtipos y tipos histológicos, entre otras.
   3.2) Análisis de grupos y subgrupos: comparar variables clínicas y de respuesta en subgrupos por subtipo, estadío o quimioterapia previa mediante tests estadísticos apropiados.
   3.3) Para ver cambios pronósticos/supervivencia se calcularán y generar curvas de supervivencia y su representación en Kaplan Meiers en toda la población y por subgrupos.
   3.4) Análisis de correlación entre las diferentes variables cuantitativas del dataset. Investigar si hay una correlación entre la edad al diagnóstico y la cantidad de mutaciones o entre el subtipo de cáncer.
   
