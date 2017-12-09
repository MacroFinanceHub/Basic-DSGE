%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico loglinealizado. 
// El Estado Estacionario es cero debido a la loglinealizacion manual.
// Modelo presentado en Clarida, Gali y Gertler (1999).
// Busqueda de malla
// (c) Carlos Rojas Quiroz 

var ygap inom pic g u;
varexo ghat uhat;
parameters alpha phi sigma psi beta lambda gamma_pic gamma_osr gamma_y 
Psi rho mu q;

phi         =0.60;
alpha       =0.50;
sigma       =1.00;
psi         =0.75; 
beta        =0.99;
rho         =0.75;
mu          =0.75;
lambda      =phi*psi*(1-(1-psi)*beta)/(1-psi);
Psi         =1/sigma;
gamma_pic   =1+(1-rho)*lambda/(rho*Psi*alpha); 
q           =1/(lambda^2+alpha*(1-beta*rho));
gamma_osr   =1.50;
gamma_y     =0.50;

model(linear);
ygap    =ygap(+1)-Psi*(inom-pic(+1))+g;
pic     =beta*pic(+1)+lambda*ygap+u;
inom    =gamma_osr*pic+gamma_y*ygap;
g       =mu*g(-1)+ghat;
u       =rho*u(-1)+uhat;
end;

options_.noprint=1;
steady;

shocks;
var ghat; stderr 1;
var uhat; stderr 1;
end;

gridd =0.05;
phipic=[0.0:gridd:2.0];
phiy  =[0.0:gridd:1.0];

for j=1:length(phipic),
    for t=1:length(phiy),
        gamma_osr  =phipic(1,j);
        gamma_y    =phiy(1,t);
        stoch_simul(irf=24,nograph, noprint);
        if info>0
                IND_PI(j,t)         = gamma_osr;
                IND_Y(j,t)          = gamma_y;
                Loss(j,t)           = 1e12;
        else 
                IND_PI(j,t)         = gamma_osr;
                IND_Y(j,t)          = gamma_y;
                Loss(j,t)           = (oo_.var(1,1)+0.5*oo_.var(3,3));
        end;
    end;
end;

Wmin=min(min(Loss));
for j=1:length(phipic),
    for t=1:length(phiy),
        if Loss(j,t)<=Wmin
                         Wmin=Loss(j,t);
                         d1min=j;
                         d2min=t;
        end;
    end;
end;

param_optimos               = [IND_PI(d1min,d2min), IND_Y(d1min,d2min)];
funcion_perdida             = Loss(d1min,d2min);
