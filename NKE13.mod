var c rnom pic lab wr y mcr a g z inv k n ce rk qr dh df da pry rer yf omega_bar s_omega;
varexo e_a e_g e_z e_yf e_s_omega;

parameters FINACC beta_C sigma h sigma_L alpha delta G_Y psi chi_inv RK Y_K epsilon I_Y K_N gamma sigma_lomega mu CE_Y C_Y nu; 
parameters theta_p chi_p rho_i phi_pic phi_y alpha_da eta etaF rho_a rho_g rho_z rho_yf rho_s_omega;
parameters RHO DEF_RATE Eomega_down ALPHA_omega GAMMA_omega a1 b1 a3 b3 a4 b4 a12 b12 a14 b14;

FINACC  = 0; //=1: under the financial accelerator; =0: w/o the financial accelerator
beta_C  = 0.99;
sigma   = 1;
h       = 0.85;
sigma_L = 1/3;
alpha   = 0.35;
delta   = 0.025;
G_Y     = 0.20;
psi     = 0.25;
chi_inv = 4.0;
RK      = ((1/beta_C)^4+0.03)^(1/4);
DEF_RATE = 0.0075;
Y_K     = (RK-(1-delta))/alpha;
epsilon = (1-delta)/RK;
I_Y     = delta/Y_K;

/**********************************************/
K_N     = 2;
B_K     = 1-1/K_N;
gamma   = 1-0.025;
[RHOSS,sigma_lomegaSS, muSS, Eomega_downSS, ALPHA_omegaSS, GAMMA_omegaSS, CE_KSS, WE_KSS, nuSS, a1SS, b1SS, a3SS, b3SS, a4SS, b4SS, a12SS, b12SS, a14SS, b14SS] = BGG_ss(beta_C,RK,delta, alpha, B_K, gamma, DEF_RATE); // Computing parameters governing the financial accelerator mechanism
RHO = RHOSS;
sigma_lomega  =sigma_lomegaSS;
mu      = muSS;
Eomega_down = Eomega_downSS;
ALPHA_omega = ALPHA_omegaSS;
GAMMA_omega = GAMMA_omegaSS;
CE_K = CE_KSS;
CE_Y      = CE_K/Y_K;
C_Y     = 1-I_Y-G_Y-CE_Y-mu*RK/Y_K*Eomega_down;
nu      = nuSS; 
a1 = a1SS;
b1 = b1SS;
a3 = a3SS;
b3 = b3SS;
a4 = a4SS;
b4 = b4SS;
a12 = a12SS;
b12 = b12SS;
a14 = a14SS;
b14 = b14SS;
/*********************************************/

theta_p = 0.75;
chi_p   = 0.95;

rho_i   = 0.90;
phi_pic = 1.5;//1.1
phi_y   = 0.125*0;

alpha_da= 0.30; 
eta     = 0.5;
etaF    = 0.5;

rho_a   = 0.99;
rho_g   = 0.95;
rho_z   = 0.00;
rho_yf  = 0.85;
rho_s_omega = 0.95;


model (linear);
// 1. Aggregate demand
da   = C_Y*c + I_Y*inv + G_Y*g + CE_Y*ce + FINACC*mu*RK/Y_K*Eomega_down*(rk+qr(-1)+k(-1) + a1*omega_bar + b1*s_omega(-1));
// 2. Aggregate demand composition - Home goods
dh   = da - eta*(pry-mcr);
// 3. Aggregate demand composition - Foreign goods
df  = da  - eta*(rer-mcr);
// 4. Euler equation of consumption 
c   = -sigma*(1-h)/(1+h)*(rnom-pic(+1)) + 1/(1+h)*c(+1) + h/(1+h)*c(-1);
// 5. Consumption of entrepreneurs
ce  = FINACC*(1-gamma)*ALPHA_omega*RK/CE_Y/Y_K*(rk+qr(-1)+k(-1) + a3*omega_bar + b3*s_omega(-1));
// 6. Spread of real return of capital over the cost of funds depend of the financial position of entrepreneurs
rk(+1) - (rnom-pic(+1)) = FINACC*(a4*omega_bar(+1) + b4*s_omega);
// 7. Real return to capital
rk  = (1-epsilon)*(pry+y-k(-1)) + epsilon*qr - qr(-1) ;
/*************************************************************************/
// 8. Tobin's Q
// qr = psi*(inv-k(-1));
// 8a. Investment adjustment cost
qr = chi_inv*(inv-inv(-1)) - beta_C*chi_inv*(inv(+1)-inv);
/*************************************************************************/
// 9. Aggregate supply
y = a + alpha*k(-1) + (1-alpha)*lab;
// 10. Supply of labor
wr  = sigma_L*lab + (1/sigma)*(1/(1-h))*c - (1/sigma)*(h/(1-h))*c(-1);
// 11. Demand for labor. Modified
wr = pry + y - lab;
// 12. Phillips curve
pic = beta_C/(1+beta_C*chi_p)*pic(+1) + chi_p/(1+beta_C*chi_p)*pic(-1) + (1-theta_p)*(1-theta_p*beta_C)/theta_p/(1+beta_C*chi_p)*mcr;
// 13. Definition of marginal cost
mcr = (1-alpha_da)*pry + alpha_da*rer;
// 14. Capital evolution
k = delta*inv + (1-delta)*k(-1);
// 15. Entrepreneur's networth evolution
n = FINACC*(gamma*ALPHA_omega*RK*K_N*(rk+qr(-1)+k(-1) + a12*omega_bar + b12*s_omega(-1)));
// 16. Monetary policy rule
rnom = rho_i*rnom(-1) + (1-rho_i)*phi_pic*pic + z;
// 16.a FIX EXCHANGE RATE
// rer = rer(-1) - pic;
// 17. Total demand for H goods
y   = (1-alpha_da)*dh + alpha_da*(yf-etaF*(pry-rer));
// 18. Perfect risk sharing
rer = (1/sigma)*(1/(1-h)*c-(h/(1-h))*c(-1)-yf); 
// 18.a Financial autarky
// pry + (yf-etaF*(pry-rer)) = rer + df;
// 19. Partcipation constraint of financial intermediaries
GAMMA_omega*RK*(FINACC*(rk+qr(-1)+k(-1)) + a14*omega_bar + FINACC*b14*s_omega(-1)) = FINACC*(1/beta_C*(qr(-1)+k(-1)) -1/beta_C/K_N*n(-1) +1/beta_C*(1-1/K_N)*(rnom(-1)-pic));
// 20. Government expenditure evolution
g   = rho_g*g(-1) + e_g;
// 21. Productivity evolution
a   = rho_a*a(-1) + e_a;
// 22. Monetary deviation evolution
z   = rho_z*z(-1) + e_z;
// 23. Foreign demand evolution
yf   = rho_yf*yf(-1) + e_yf;
// 24. risk shock
s_omega = rho_s_omega*s_omega(-1) + e_s_omega;
end;

initval;
c       = 0;
rnom    = 0;
pic     = 0;
g       = 0;
wr      = 0;
lab     = 0;
y       = 0;
a       = 0;
mcr     = 0;
z       = 0;
inv     = 0;
k       = 0;
n       = 0;
ce      = 0;
rk      = 0;
qr      = 0;
dh      = 0;
df      = 0;
da      = 0; 
pry     = 0;
rer     = 0;
yf      = 0;
omega_bar = 0;
s_omega   = 0;
end;

resid;

shocks;
var e_a; stderr 0.7;
var e_g; stderr 0.5;
var e_z; stderr 0.25;
var e_yf; stderr 2.0;
var e_s_omega; stderr 8.0;
end;

steady;

check;

stoch_simul(order=1,irf=24,nograph);