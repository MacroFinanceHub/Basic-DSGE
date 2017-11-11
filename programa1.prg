'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
' Programa para obtener estadísticos de los componentes cíclicos de las siguientes cuentas:
' PBI real, Empleo (10 a más trabajadores), Gasto Público (consumo más inversión pública) e Inversión privada
'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wfcreate(page=Ciclos) LambdaDSGE q 1980.1 2016.4
import "C:\Users\crojas\Documents\00.Personales\Lambda\Analisis Macro con DSGEs\Cap. 1\BaseLambda01.xlsx" "Data" @freq q 1980.1

smpl @all
' Desestacionalización de las variables de interés
yreal.x12
lab.x12
gpub.x12
inv.x12
kap.x12
cpriv.x12

series yreal=yreal_sa
series lab=lab_sa
series gpub=gpub_sa
series inv=inv_sa
series kap=kap_sa
series cpriv=cpriv_sa
delete yreal_sa lab_sa gpub_sa inv_sa kap_sa cpriv_sa


for %x yreal lab kap gpub inv cpriv
' Prediciendo doce períodos hacia adelante para evitar problema de colas
range @first 2019.4
equation eqf_{%x}.ls {%x} c ar(1)
smpl 2017.1 2019.4
eqf_{%x}.forecast {%x}_f
' Obteniendo componentes cíclicos
smpl 1999.1 2019.4
series log_{%x}=log({%x}_f)
' Hodrick-Prescott
log_{%x}.hpf hp_trend_{%x}
series hp_cyc_{%x}= log_{%x} - hp_trend_{%x}
delete hp_trend_{%x}
' Baxter y King
log_{%x}.bpf bk_cyc_{%x}
' Tendencia cuadrática
equation eq_{%x}.ls log({%x}) c @trend @trend^2
eq_{%x}.makeresids tc_cyc_{%x}
range @first 2016.4
next

' Obteniendo la PTF cíclica a partir del residuo de solow
scalar alpha=0.67
smpl 2002.1 2016.4
for %y hp bk tc
series {%y}_ptf={%y}_cyc_yreal - alpha*{%y}_cyc_lab - (1-alpha)*{%y}_cyc_kap
next 

' Construyendo tablas de resultados
for %m tc bk hp

' Tabla 1: Diversos estadísticos bajo diferentes filtros 
table stats_{%m}

' Variables bajo consideración
stats_{%m}(2,1)="PBI"
stats_{%m}(3,1)="Empleo"
stats_{%m}(4,1)="Inversion"
stats_{%m}(5,1)="Gasto publico"
stats_{%m}(6,1)="PTF"
stats_{%m}(7,1)="Consumo"

' Estadísticos bajo consideración
stats_{%m}(1,2)="D.E"
stats_{%m}(1,3)="D.E. relativa PBI"
stats_{%m}(1,4)="Corr. PBI"
stats_{%m}(1,5)="Autocorr"

' Cálculo correspondiente
stats_{%m}(2,2)=@stdev({%m}_cyc_yreal)
stats_{%m}(2,3)=@stdev({%m}_cyc_yreal)/@stdev({%m}_cyc_yreal)
stats_{%m}(2,4)=@cor({%m}_cyc_yreal,{%m}_cyc_yreal)
stats_{%m}(2,5)=@cor({%m}_cyc_yreal,{%m}_cyc_yreal(-1))

stats_{%m}(3,2)=@stdev({%m}_cyc_lab)
stats_{%m}(3,3)=@stdev({%m}_cyc_lab)/@stdev({%m}_cyc_yreal)
stats_{%m}(3,4)=@cor({%m}_cyc_lab,{%m}_cyc_yreal)
stats_{%m}(3,5)=@cor({%m}_cyc_lab,{%m}_cyc_lab(-1))

stats_{%m}(4,2)=@stdev({%m}_cyc_inv)
stats_{%m}(4,3)=@stdev({%m}_cyc_inv)/@stdev({%m}_cyc_yreal)
stats_{%m}(4,4)=@cor({%m}_cyc_inv,{%m}_cyc_yreal)
stats_{%m}(4,5)=@cor({%m}_cyc_inv,{%m}_cyc_inv(-1))

stats_{%m}(5,2)=@stdev({%m}_cyc_gpub)
stats_{%m}(5,3)=@stdev({%m}_cyc_gpub)/@stdev({%m}_cyc_yreal)
stats_{%m}(5,4)=@cor({%m}_cyc_gpub,{%m}_cyc_yreal)
stats_{%m}(5,5)=@cor({%m}_cyc_gpub,{%m}_cyc_gpub(-1))

stats_{%m}(6,2)=@stdev({%m}_ptf)
stats_{%m}(6,3)=@stdev({%m}_ptf)/@stdev({%m}_cyc_yreal)
stats_{%m}(6,4)=@cor({%m}_ptf,{%m}_cyc_yreal)
stats_{%m}(6,5)=@cor({%m}_ptf,{%m}_ptf(-1))

stats_{%m}(7,2)=@stdev({%m}_cyc_cpriv)
stats_{%m}(7,3)=@stdev({%m}_cyc_cpriv)/@stdev({%m}_cyc_yreal)
stats_{%m}(7,4)=@cor({%m}_cyc_cpriv,{%m}_cyc_yreal)
stats_{%m}(7,5)=@cor({%m}_cyc_cpriv,{%m}_cyc_cpriv(-1))
next

' Tabla 2: Parámetros de procesos exógenos bajo diferentes filtros 
table param

param(2,1)="rho_z"
param(3,1)="sigma_z"
param(4,1)="rho_g"
param(5,1)="sigma_g"
param(1,2)="Hodrick-Prescott"
param(1,3)="Baxter-King"
param(1,4)="Tendencia cuadrática"

param(2,2)=@cov(hp_ptf,hp_ptf(-1))/(@stdev(hp_ptf)^2)
param(3,2)=((1-(@cov(hp_ptf,hp_ptf(-1))/(@stdev(hp_ptf)^2))^2)*(@stdev(hp_ptf)^2))^0.5
param(4,2)=@cov(hp_cyc_gpub,hp_cyc_gpub(-1))/(@stdev(hp_cyc_gpub)^2)
param(5,2)=((1-(@cov(hp_cyc_gpub,hp_cyc_gpub(-1))/(@stdev(hp_cyc_gpub)^2))^2)*(@stdev(hp_cyc_gpub)^2))^0.5

param(2,3)=@cov(bk_ptf,bk_ptf(-1))/(@stdev(bk_ptf)^2)
param(3,3)=((1-(@cov(bk_ptf,bk_ptf(-1))/(@stdev(bk_ptf)^2))^2)*(@stdev(bk_ptf)^2))^0.5
param(4,3)=@cov(bk_cyc_gpub,bk_cyc_gpub(-1))/(@stdev(bk_cyc_gpub)^2)
param(5,3)=((1-(@cov(bk_cyc_gpub,bk_cyc_gpub(-1))/(@stdev(bk_cyc_gpub)^2))^2)*(@stdev(bk_cyc_gpub)^2))^0.5

param(2,4)=@cov(tc_ptf,tc_ptf(-1))/(@stdev(tc_ptf)^2)
param(3,4)=((1-(@cov(tc_ptf,tc_ptf(-1))/(@stdev(tc_ptf)^2))^2)*(@stdev(tc_ptf)^2))^0.5
param(4,4)=@cov(tc_cyc_gpub,tc_cyc_gpub(-1))/(@stdev(tc_cyc_gpub)^2)
param(5,4)=((1-(@cov(tc_cyc_gpub,tc_cyc_gpub(-1))/(@stdev(tc_cyc_gpub)^2))^2)*(@stdev(tc_cyc_gpub)^2))^0.5


