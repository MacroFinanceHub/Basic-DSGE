%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - NKE %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico loglinealizado. 
// Economia abierta.
// Fixed exchange rate
// (c) Carlos Rojas Quiroz 

var c ch cf prh rer rnom pic pih lab wr y mchr a g z yf dep;
varexo e_a e_g e_z e_yf;

parameters beta_C sigma alpha_c eta_c sigma_L theta_p phi_pic phi_y eta_f rho_a rho_g rho_z rho_yf;

beta_C  = 0.995;
sigma   = 1;
eta_c   = 1;
eta_f   = 1;
alpha_c = 0.30;
sigma_L = 1.0;
theta_p = 0.75;
phi_pic = 1.5;
phi_y   = 0.125;

rho_a   = 0.90;
rho_g   = 0.80;
rho_z   = 0.70;
rho_yf  = 0.85;


model (linear);
ch  = c - eta_c*prh;
cf  = c - eta_c*rer;
0   = (1-alpha_c)*prh + alpha_c*rer;
c   = -sigma*(rnom-pic(+1)) + c(+1) + sigma*g*(1-rho_g);
wr  = sigma_L*lab + (1/sigma)*c;
y   = a + lab;
pih = beta_C*pih(+1) + (1-theta_p)*(1-theta_p*beta_C)/theta_p*(mchr-prh);
mchr= wr - a;
dep = 0; 
y   = (1-alpha_c)*ch + alpha_c*(yf-eta_f*(prh-rer));
pic = pih + alpha_c/(1-alpha_c)*(rer-rer(-1));
rer = (1/sigma)*(c-yf) -g;
a   = rho_a*a(-1) + e_a;
g   = rho_g*g(-1) + e_g;
z   = rho_z*z(-1) + e_z;
yf  = rho_yf*yf(-1) + e_yf;
dep = rer-rer(-1) + pic;
end;

initval;
c       = 0;
ch      = 0;
cf      = 0;
rnom    = 0;
pic     = 0;
pih     = 0;
rer     = 0;
prh     = 0;
g       = 0;
wr      = 0;
lab     = 0;
y       = 0;
a       = 0;
mchr    = 0;
z       = 0;
yf      = 0;
dep     = 0;
end;

resid;

shocks;
var e_a; stderr 0.7;
var e_g; stderr 0.5;
var e_z; stderr 0.25;
var e_yf; stderr 1.0;
end;

steady;

check;

stoch_simul(order=1,irf=24,nograph);