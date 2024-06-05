# cbioportal-project-UB
Repositorio para el trabajo final del curso de Ciencia de Datos de IL-3 Universidad de Barcelona.

Pasos:
1) Buscar un dataset en cBioportal con las características que nos interesan. https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018
2) Conectar cBioportal con R studio mediante una interfaz-API para obtener el dataset escogido.
3) Explotar los datos con R studio:
   3.1) Generar tablas con estadísticas descriptivas para la mayoría de variables con el cálculo de medidas de tendencia central para variables como: edad al diagnóstico, alteraciones genómicas, subtipos y tipos histológicos, entre otras; utilizando paquetes como “gtsummary” y la funnción tbl_summary().
   3.2) Análisis de grupos y subgrupos: comparar variables clínicas y de respuesta en subgrupos por subtipo, estadío o quimioterapia previa mediante tests estadísticos apropiados.
   3.3) Para ver cambios pronósticos/supervivencia se calcularán y generar curvas de supervivencia y su representación en Kaplan Meiers en toda la población y por subgrupos, mediante los paquetes: “survival” y “survminer” ejecutando funciones como survfit o ggsurvplot(), respectivamente.
   3.4) Análisis de correlación entre las diferentes variables cuantitativas del dataset. Investigar si hay una correlación entre la edad al diagnóstico y la cantidad de mutaciones o entre el subtipo de cáncer y la expresión de mRNA**.
5) Crear alguna interfaz con Shiny app para visualizar de forma interactiva los gráficos de supervivencia KM por subgrupos de variables pronósticas que encontremos en el análisis exploratorio así como tablas de prevalencia de mutaciones en función del subtipo escogido, estadío o tipo histológico.

Git y Github en RStudio: https://rpubs.com/RonaldoAnticona/818156
