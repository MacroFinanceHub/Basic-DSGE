# Welcome!
This site contains some basics on DSGE modeling. It's part of workshops and classes I teach for undergraduate students (National Engineering University and Peruvian University of Applied Sciences). The platform I use is Dynare/Matlab.

# Real business cycles 
All the next models are composed of three agents: consumers, producers, and government (with a simple exogenous process). The calibration follows Fernandez-Villaverde's notes on Macroeconomics. 

RBC01.mod:       Nonlinear First Order Conditions (FOC). Linearized with a first order Taylor expansion.

RBC01b.mod:      Nonlinear FOC, log-linearized using a first order Taylor expansion.

RBC01b_graph.mod:RBC01b with changes in the productivity's autoregressive parameter. 

RBC02.mod:       Hand-operated log-linearization of equilibrium conditions. 

RBC02b.mod:      Model with a fixed supply of labor (substitution and income effects cancel each other).

RBC03.mod:       Greenwood-Hercowitz-Huffman utility function.

RBC04.mod:       King-Plosser-Rebelo utility function.

RBC05.mod:       Constant Risk Relative Aversion coefficient utility function.

RBC06.mod:       Internal consumption habits.

RBC07.mod:       External consumption habits.

RBC08.mod:       Capital utilization rate as a variable cost.

RBC09.mod:       Capital utilization rate as a depreciation rate.

RBC10.mod:       Investment's adjusment cost. 

RBC11.mod:       Capital's adjustment cost. 

RBC12.mod:       Specific investment shock. 

RBC13.mod:       Hansen's model with indivisible Labor.  

RBC14.mod:       Cash-in-Advance constraint.

RBC15.mod:       Money in Utility.

RBC16.mod:       Small open economy - variable internal discount factor.

RBC16a.mod:      Small open economy - variable external discount factor.

RBC17.mod:       Small open economy - debt elastic interest rate.

RBC18.mod:       Small open economy - portfolio's adjustment cost.

RBC19.mod:       Small open economy - complete markets.

RBC20.mod:       Markup prices shock.

RBC21_Can.mod:   Small open economy - trend shock - Canada's calibration.

RBC21_Mex.mod:   Small open economy - trend shock - Mexico's calibration.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Grafico_program.m:      Program for Impulse-Response graphics (IRF).

Grafico_programSOE.m :  Replicating Figure 1 from SG-U (2003).

Grafico_programSOE2.m:  Replicating Figure 3 from A-G (2007).

Grid_alpha.md:          Setting different values of alpha under a productivity shock.

BaseLambda01.xlsx:      Base de datos para obtención de ciclos económicos.

programa1.prg:   Programa en Eviews para obtención de ciclos económicos.

RBCTAREA.mod:    Programa que incorpora choque de preferencias y de oferta laboral al modelo RBC básico.

TAREA_graph.m:   Programa para graficar los dos choques de la Tarea 1.

NKE01.mod:       Modelo neokeynesiano de tres ecuaciones con politica monetaria optima (manual) extraido de C-G-G (1999).

NKE02.mod:       Modelo neokeynesiano de tres ecuaciones con politica monetaria optima - comando osr.

NKE03.mod:       Modelo neokeynesiano de tres ecuaciones con politica monetaria optima - comando discretionary_policy.

NKE04.mod:       Modelo neokeynesiano de tres ecuaciones con politica monetaria optima - comando ramsey_policy.

NKE05.mod:       Modelo neokeynesiano de tres ecuaciones con politica monetaria optima - frontera de politica.

NKE06.mod:       Modelo neokeynesiano semiestructural backward-looking.

NKE07.mod:       Modelo neokeynesiano en niveles.

NKE08.mod:       Modelo neokeynesiano log-linealizado.

NKE09.mod:       Modelo neokeynesiano log-linealizado incluyendo inversion y gasto publico. 

NKE10.mod:       Modelo neokeynesiano estándar para economia pequena y abierta.

NKE11.mod:       Modelo neokeynesiano de tres ecuaciones con politica monetaria optima - busqueda de malla.

GM01.mod:        Economia pequena y abierta - regla sobre inflacion domestica.

GM02.mod:        Economia pequena y abierta - regla sobre inflacion IPC.

GM03.mod:        Economia pequena y abierta - tipo de cambio nominal fijo.

GM04.mod:        Economia pequena y abierta - regla optima.

Graph_NKE.mod:   Programa para replicar la Figura 1 de G-M (2005).
