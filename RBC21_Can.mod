%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico para economia pequena y abierta
// Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// Modelo calibrado para Canada con trend shock.
// (c) Carlos Rojas Quiroz 

var c k y b q g l u z uc ul c_y i_y invest nx i_y_percentage c_y_percentage 
    log_y log_c log_i
    delta_y;

predetermined_variables k b;

varexo eps_z eps_g;

parameters mu_g sigma rho_g delta phi psi b_star alpha rho_z r_star beta gamma b_share;

% Valores de la Tabla 3
beta    = 1/1.02;
gamma   = 0.36;
b_share = 0.1; 
psi     = 0.001;
alpha   = 0.68;
sigma   = 2;
delta   = 0.05;
phi     = 4;

% Parametros estimados (Tabla 4)
mu_g = log(1.0073);
rho_z = 0.95;
rho_g = 0.01;
 
model;
% Funcion de Produccion
y = exp(z)*k^(1-alpha)*(exp(g)*l)^alpha; 
% Choque transitorio 
z = rho_z*z(-1)+eps_z; 
% Choque permanente 
g = (1-rho_g)*mu_g+rho_g*g(-1)+eps_g; 
% Funcion de utilidad 
u = (c^gamma*(1-l)^(1-gamma))^(1-sigma)/(1-sigma); 
% Utilidad marginal del consumo
uc = gamma*u/c*(1-sigma); 
% Desutilidad marginal del trabajo
ul = -(1-gamma)*u/(1-l)*(1-sigma); 
% Restriccion presupuestaria
c+exp(g)*k(+1)=y+(1-delta)*k-phi/2*(exp(g)*k(+1)/k-exp(mu_g))^2*k-b+q*exp(g)*b(+1); 
% Precio de la deuda
1/q = 1+r_star+psi*(exp(b(+1)-b_star)-1);
% Condicion de primer orden respecto al Capital
uc*(1+phi*(exp(g)*k(+1)/k-exp(mu_g)))*exp(g)=beta*exp(g*(gamma*(1-sigma)))*uc(+1)*(1-delta+(1-alpha)*y(+1)/k(+1)
                                                  -phi/2*(2*(exp(g(+1))*k(+2)/k(+1)-exp(mu_g))*(-1)*exp(g(+1))*k(+2)/k(+1)+
                                                          (exp(g(+1))*k(+2)/k(+1)-exp(mu_g))^2));
% Condicion de primer orden respecto al trabajo
ul+uc*alpha*y/l=0; 
% Ecuacion de Euler para los bonos
uc*exp(g)*q=beta*exp(g*(gamma*(1-sigma)))*uc(+1); 
% Ecuacion de movimiento del capital
invest = exp(g)*k(+1)-(1-delta)*k+phi/2*(exp(g)*k(+1)/k-exp(mu_g))^2*k; 

% Definiciones de variables para los IRFs

% Ratio Consumo/PBI
c_y = c/y;
% Ratio Inversion/PBI
i_y = invest/y;
% Ratio Exportaciones Netas/PBI
nx=(b-exp(g)*q*b(+1))/y;

% Variables auxiliares
i_y_percentage=log(i_y);
c_y_percentage=log(c_y);
log_y=log(y);
log_c=log(c);
log_i=log(invest);

% Definicion de tasa de crecimiento del PBI
delta_y=log(y)-log(y(-1))+g(-1);

end;

steady_state_model;
q=beta*exp(mu_g)^(gamma*(1-sigma)-1);
YKbar=((1/q)-(1-delta))/(1-alpha);
c_y=1+(1-exp(mu_g)-delta)*(1/YKbar)-(1-exp(mu_g)*q)*b_share;
l=(alpha*gamma)/(c_y-gamma*c_y+alpha*gamma);
k=(((exp(mu_g)^alpha)*(l^alpha))/YKbar)^(1/alpha);
y=k^(1-alpha)*(l*exp(mu_g))^alpha;
c=c_y*y;
invest=(exp(mu_g)-1+delta)*k;
b_star=b_share*y;
nx=(y-c-invest)/y;
r_star = 1/q-1;

b = b_star;
z = 0;
g = mu_g;
u = (c^gamma*(1-l)^(1-gamma))^(1-sigma)/(1-sigma);
uc = gamma*u/c*(1-sigma);
ul = -(1-gamma)*u/(1-l)*(1-sigma);
i_y = (exp(g)*k-(1-delta)*k)/y;
i_y_percentage=log(i_y);
c_y_percentage=log(c_y);
log_y=log(y);
log_c=log(c);
log_i=log(invest);
delta_y=mu_g;
end;

shocks;
var eps_g; stderr 0.01; 
var eps_z; stderr 0.01;
end;

resid;
steady;
check;

stoch_simul(order=1,irf=26,nograph);
