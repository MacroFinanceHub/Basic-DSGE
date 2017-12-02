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

RBC16.mod:       Economia pequena y abierta - factor de descuento endogeno.

RBC16a.mod:      Economia pequena y abierta - factor de descuento endogeno sin internalizacion.

RBC17.mod:       Economia pequena y abierta - tasa de interes elastica a la deuda.

RBC18.mod:       Economia pequena y abierta - costo de ajuste del portafolio.

RBC19.mod:       Economia pequena y abierta - mercados completos.

RBC20.mod:       Función de Utilidad logarítmica con choque al margen de precios.

RBC21_Can.mod:   Economia pequena y abierta - trend shock - calibracion para Canada

RBC21_Mex.mod:   Economia pequena y abierta - trend shock - calibracion para Mexico.

Grafico_program.m: Programa para hacer gráficos en Matlab a partir de productos del Dynare (IRF).

Grafico_programSOE.m : Programa para replicar la Figura 1 de SG-U (2003).

Grafico_programSOE2.m: Programa para replicar la Figura 3 de A-G (2007).

Grid_alpha.md:   "Setea" valores de alpha ante choques de productividad.

BaseLambda01.xlsx: Base de datos para obtención de ciclos económicos.

programa1.prg:   Programa en Eviews para obtención de ciclos económicos.

RBCTAREA.mod:    Programa que incorpora choque de preferencias y de oferta laboral al modelo RBC básico.

TAREA_graph.m:   Programa para graficar los dos choques de la Tarea 1.
