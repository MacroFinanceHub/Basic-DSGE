# Topicos-DSGE
Modelos RBC y NK en Dynare 
Los modelos utilizan el entorno Dynare para Matlab. Todos ellos estan compuestos de Consumidores, Productores y Gobierno que actuan en competencia perfecta. Los consumidores son dueños del capital. La calibración sigue a Fernandez-Villaverde (Lecture Notes on Macroeconomics). 

RBC01.mod:       Condiciones de primer orden (CPO) no lineales, linealizadas mediante una expansion de Taylor de primer orden. 
RBC01b.mod:      CPO no lineales, log-linealizadas mediante una expansion de Taylor de primer orden. 
RBC01b_graph.mod:RBC01b con cambios en el parámetro de persistencia del choque de productividad.
RBC02.mod:       CPO log-linealizadas manualmente. 
RBC02b.mod:      Oferta de trabajo fija (efectos sustitución e ingreso se cancelan).
RBC03.mod:       Función de Utilidad a la Greenwood-Hercowitz-Huffman.
RBC04.mod:       Función de Utilidad a la King-Plosser-Rebelo.
RBC05.mod:       Función de Utilidad con Coeficiente de Aversión Relativa al Riesgo constante.
RBC06.mod:       Función de Utilidad logarítmica con hábitos de consumo interno.
RBC07.mod:       Función de Utilidad logarítmica con hábitos de consumo externo.
RBC08.mod:       Función de Utilidad logarítmica con ratio de uso de capital en forma de costo.
RBC09.mod:       Función de Utilidad logarítmica con ratio de uso de capital con depreciación variable. 
RBC10.mod:       Función de Utilidad logarítmica con costos de ajuste a la inversión. 
RBC11.mod:       Función de Utilidad logarítmica con costos de ajuste al capital. 
RBC12.mod:       Función de Utilidad logarítmica con choque específico a la inversión. 
RBC13.mod:       Función de Utilidad logarítmica con dinero en una restricción Cash-in-Advance.  
RBC14.mod:       Función de Utilidad logarítmica con dinero en la función de utilidad.
RBC15.mod:       Función de Utilidad logarítmica con trabajo indivisible a la Hansen.
RBC20.mod:       Función de Utilidad logarítmica con choque al margen de precios.
Grafico_program.m: Programa para hacer gráficos en Matlab a partir de productos del Dynare (IRF)
Grid_alpha.md:   "Setea" valores de alpha ante choques de productividad.
