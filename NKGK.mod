%----------------------------------------------------------------
% Modelo BASE: BGG + GK (sin regla macroprudencial)
%----------------------------------------------------------------
 
close all;

%----------------------------------------------------------------
% 1. Preámbulo
%----------------------------------------------------------------

// Variables endógenas
var Welfare Welfare_bank lam c h hd w wtil mcW fW DelW inv k rK q y yC yF yH xF xH xHstar R xi RL pic rer pH ptilH pF ptilF pY piS mcH fH DelH mcF fF DelF dstar m tb stb sdstar gam_Y gam_L gam_C gam_I gam_W gam_K gam_YC gam_YF gam_YH gam_XF gam_XH gam_XHstar gam_M gam_LAM piH piF
varrhoL varrhoN LEV L d n Util;

// Variables de estado exógenas
var v u z a zeta eR yCo Rstar pistar pCostar ystar g fin;
   
// Innovaciones exógenas
varexo eps_v eps_u eps_z eps_a eps_zeta eps_eR eps_yCo eps_Rstar eps_pistar eps_pCostar eps_ystar eps_g eps_fin;

// Parámetros    
parameters SIGMA PHI ALPHA DELTA EPSILON_H EPSILON_F EPSILON_W O CHI;
parameters VARSIGMA PSI ETA ETASTAR GAMA THETA_W VARTHETA_W THETA_H VARTHETA_H THETA_F VARTHETA_F RHO_R ALPHA_PI ALPHA_Y ALPHA_F;
parameters BETTA OSTAR KAPPA;
parameters RHO_v RHO_u RHO_z RHO_a RHO_zeta RHO_eR RHO_yCo RHO_Rstar RHO_pistar RHO_pCostar RHO_ystar RHO_g RHO_fin;
parameters SIG_v SIG_u SIG_z SIG_a SIG_zeta SIG_eR SIG_yCo SIG_Rstar SIG_pistar SIG_pCostar SIG_ystar SIG_g SIG_fin;
parameters lam_ss c_ss h_ss hd_ss w_ss wtil_ss mcW_ss fW_ss DelW_ss inv_ss k_ss rK_ss q_ss y_ss yC_ss yF_ss yH_ss xF_ss xH_ss xHstar_ss R_ss xi_ss RL_ss pi_ss rer_ss piS_ss pH_ss ptilH_ss pF_ss ptilF_ss pY_ss mcH_ss fH_ss DelH_ss mcF_ss fF_ss DelF_ss dstar_ss m_ss tb_ss stb_ss sdstar_ss;
parameters v_ss u_ss z_ss a_ss zeta_ss eR_ss yCo_ss Rstar_ss pistar_ss pCostar_ss ystar_ss g_ss fin_ss;
parameters sCo_ss sg_ss;
parameters GAMMA_ss LEV_ss iota_ss alphaLK ww_ss varrhoN_ss varrhoL_ss mu_ss n_ss d_ss Util_ss Welfare_ss Welfare_bank_ss;

// Parámetros calibrados (Modelo neokeynesiano simple)
SIGMA=1;             // log utility (Medina and Soto, 2007)
PHI=1;               // unitary Frisch elasticity (Adolfson et al., 2008)
ALPHA=1-0.66;        // labor share of 66% (Medina and Soto, 2007)
DELTA=0.06/4;        // annual depreciation rate of 6% (Medina and Soto, 2007)
EPSILON_H=11;        // steady state price markup of 10% (Medina and Soto, 2007)
EPSILON_F=11;        // steady state price markup of 10% (Medina and Soto, 2007)
EPSILON_W=11;        // steady state wage markup of 10% (Medina and Soto, 2007)
O=0.32;              // home bias in domestic demand of 68% (imports/domestic demand=32%, 1987-2012)
CHI=0.61;            // CHI=c+(1-c)*t, c=0.4 (production Codelco/total, 1987-2012), t = 0.35 (general tax)

// Parámetros calibrados (Empresarios y Bancos)
GAMMA_ss=1.0380^0.25-1;
LEV_ss=13;
iota_ss=0.002;
alphaLK=0.51;

// Valores de Estado Estacionario (para Chile)
stb_ss=0.04;         // average share of trade balance / GDP, 1987-2012
sg_ss=0.11;          // average share of government consumption / GDP, 1987-2012
sCo_ss=0.10;         // average share of GDP copper sector / total, 1987-2012
pi_ss=1.03^.25;      // central bank inflation target, 2001-2012
a_ss=1.025^.25;      // quarterly balanced growth path (average growth rate 4.5%-2.0% approx. labor force growth, 2001-2012)
R_ss=1.058^.25;      // quarterly gross home nominal interest rate (Fuentes and Gredig, 2008)
Rstar_ss=1.045^.25;  // quarterly gross foreign nominal interest rate (Fuentes and Gredig, 2008)
xi_ss=1.014^.25;     // average quarterly gross country (EMBI Chile) spread, 2001-2012
h_ss=0.2;            // normalization
v_ss=1;              // normalization
u_ss=1;              // normalization
pH_ss=1;             // normalization
z_ss=1;              // normalization 
zeta_ss=1;           // normalization
eR_ss=1;             // normalization
ystar_ss=1;          // normalization
pCostar_ss=1;        // normalization
fin_ss=1;            // normalization

// Parámetros estimados (García-Cicco et al. 2014)
VARSIGMA=0.70;
PSI=0.008;
ETA=1.30; 
ETASTAR=0.32;
GAMA=0.43;
THETA_W=0.94;
VARTHETA_W=0.37; 
THETA_H=0.46;
VARTHETA_H=0.39;
THETA_F=0.81;
VARTHETA_F=0.39;

// Parámetros de política (regla monetaria) 
RHO_R   =0.750; 
ALPHA_PI=1.500; 
ALPHA_Y =0.110;
ALPHA_F =0.000;

// Parámetros de persistencia de choques exógenos
RHO_v=0.76;
RHO_u=0.74; 
RHO_z=0.74; 
RHO_a=0.35;
RHO_zeta=0.92;
RHO_eR=0.0;
RHO_yCo=0.4794; 
RHO_g=0.6973;
RHO_Rstar=0.9614; 
RHO_ystar=0.8665; 
RHO_pistar=0.3643; 
RHO_pCostar=0.9620;
RHO_fin=0.85;

// Magnitud de los choques
SIG_v=0.022;  
SIG_u=0.021; 
SIG_z=0.010;
SIG_a=0.003; 
SIG_zeta=0.001;
SIG_eR=0.002;
SIG_yCo=0.0293; 
SIG_g=0.0145; 
SIG_Rstar=0.0011; 
SIG_ystar=0.0062; 
SIG_pistar=0.0273;
SIG_pCostar=0.1413;
SIG_fin=0.01;

%----------------------------------------------------------------
% 2. Modelo (85=38+13+18+13+3 ecuaciones)
%----------------------------------------------------------------

model; 

// Condiciones de equilibrio: Modelo Neokeynesiano Simple (38)

exp(lam)=(exp(c)-VARSIGMA*exp(c(-1))/exp(a(-1)))^(-SIGMA)-BETTA*VARSIGMA*exp(v(+1))/exp(v)*(exp(c(+1))*exp(a)-VARSIGMA*exp(c))^(-SIGMA); //E1
exp(w)*exp(mcW)=KAPPA*exp(h)^PHI/exp(lam); //E2
exp(lam)=BETTA/exp(a)*exp(R)*exp(v(+1))/exp(v)*exp(lam(+1))/exp(pic(+1)); //E3
exp(lam)=BETTA/exp(a)*exp(Rstar)*exp(xi)*exp(v(+1))/exp(v)*exp(piS(+1))*exp(lam(+1))/exp(pic(+1)); //E4
exp(yC)=((1-O)^(1/ETA)*exp(xH)^((ETA-1)/ETA)+O^(1/ETA)*exp(xF)^((ETA-1)/ETA))^(ETA/(ETA-1)); //E5
exp(xF)=O*exp(pF)^(-ETA)*exp(yC); //E6
exp(xH)=(1-O)*exp(pH)^(-ETA)*exp(yC); //E7
exp(mcH)=exp(rK)^ALPHA*exp(w)^(1-ALPHA)/(exp(pH)*exp(z)*exp(a)^(1-ALPHA))/(ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)); //E8
exp(fH)=exp(ptilH)^-EPSILON_H*exp(yH)*exp(mcH)+BETTA*THETA_H*exp(v(+1))/exp(v)*(exp(lam(+1))/exp(lam)*(exp(pic)^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pic(+1)))^-EPSILON_H*(exp(ptilH)/exp(ptilH(+1)))^-EPSILON_H*(exp(pH)/exp(pH(+1)))^(-1-EPSILON_H)*exp(fH(+1))); //E9
exp(fH)=exp(ptilH)^(1-EPSILON_H)*exp(yH)*(EPSILON_H-1)/EPSILON_H+BETTA*THETA_H*exp(v(+1))/exp(v)*(exp(lam(+1))/exp(lam)*(exp(pic)^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pic(+1)))^(1-EPSILON_H)*(exp(ptilH)/exp(ptilH(+1)))^(1-EPSILON_H)*(exp(pH)/exp(pH(+1)))^-EPSILON_H*exp(fH(+1))); //E10
exp(xHstar)=OSTAR*(exp(pH)/exp(rer))^-ETASTAR*exp(ystar); //E11
exp(R)/R_ss=(exp(R(-1))/R_ss)^RHO_R*((exp(pic)/pi_ss)^ALPHA_PI*(exp(y)/exp(y(-1)))^ALPHA_Y*(exp(L)/exp(L(-1)))^ALPHA_F)^(1-RHO_R)*exp(eR); //E12
exp(yH)*exp(DelH)=exp(z)*(exp(k(-1))/exp(a(-1)))^ALPHA*(exp(a)*exp(hd))^(1-ALPHA); //E13
1=(1-THETA_H)*exp(ptilH)^(1-EPSILON_H)+THETA_H*(exp(pH(-1))/exp(pH)*exp(pic(-1))^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pic))^(1-EPSILON_H); //E14
exp(yH)=exp(xH)+exp(xHstar); //E15
exp(yC)=exp(c)+exp(inv)+exp(g); //E16
exp(rer)/exp(rer(-1))=exp(piS)*exp(pistar)/exp(pic); //E17
exp(y)=exp(c)+exp(inv)+exp(g)+exp(xHstar)+exp(yCo)-exp(m); //E18
tb=exp(pH)*exp(xHstar)-exp(rer)*exp(m)+exp(rer)*exp(pCostar)*exp(yCo); //E19
exp(rer)*dstar=exp(rer)*dstar(-1)/exp(a(-1))/exp(pistar)*exp(Rstar(-1))*exp(xi(-1))-tb+(1-CHI)*exp(rer)*exp(pCostar)*exp(yCo); //E20
exp(xi)=xi_ss*exp(PSI*(exp(rer)*dstar-rer_ss*dstar_ss)/(rer_ss*dstar_ss)+(exp(zeta)-zeta_ss)/zeta_ss); //E21, for calibrations with dstar_ss<0, write exp(PSI*(exp(rer)*dstar-rer_ss*dstar_ss)/(-rer_ss*dstar_ss)...)
exp(DelH)=(1-THETA_H)*exp(ptilH)^-EPSILON_H+THETA_H*(exp(pH(-1))/exp(pH)*exp(pic(-1))^VARTHETA_H*pi_ss^(1-VARTHETA_H)/exp(pic))^-EPSILON_H*exp(DelH(-1)); //E22
exp(k)=(1-DELTA)*exp(k(-1))/exp(a(-1))+(1-GAMA/2*(exp(inv)/exp(inv(-1))*exp(a(-1))-a_ss)^2)*exp(u)*exp(inv); //E23
exp(k(-1))/exp(hd)=exp(a(-1))*ALPHA/(1-ALPHA)*exp(w)/exp(rK); //E24
1/exp(q)=(1-GAMA/2*(exp(inv)/exp(inv(-1))*exp(a(-1))-a_ss)^2-GAMA*(exp(inv)/exp(inv(-1))*exp(a(-1))-a_ss)*exp(inv)/exp(inv(-1))*exp(a(-1)))*exp(u)+BETTA/exp(a)*GAMA*exp(v(+1))/exp(v)*exp(q(+1))/exp(q)*exp(lam(+1))/exp(lam)*(exp(inv(+1))/exp(inv)*exp(a)-a_ss)*(exp(inv(+1))/exp(inv)*exp(a))^2*exp(u(+1)); //E25
exp(pY)*exp(y)=exp(c)+exp(inv)+exp(g)+tb; //E26
1=(1-THETA_F)*exp(ptilF)^(1-EPSILON_F)+THETA_F*(exp(pF(-1))/exp(pF)*exp(pic(-1))^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pic))^(1-EPSILON_F); //E27
exp(yF)=exp(xF); //E28
exp(m)=exp(yF)*exp(DelF); //E29
exp(DelF)=(1-THETA_F)*exp(ptilF)^-EPSILON_F+THETA_F*(exp(pF(-1))/exp(pF)*exp(pic(-1))^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pic))^-EPSILON_F*exp(DelF(-1)); //E30
exp(mcF)=exp(rer)/exp(pF); //E31
exp(fF)=exp(ptilF)^-EPSILON_F*exp(yF)*exp(mcF)+BETTA*THETA_F*exp(v(+1))/exp(v)*(exp(lam(+1))/exp(lam)*(exp(pic)^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pic(+1)))^-EPSILON_F*(exp(ptilF)/exp(ptilF(+1)))^-EPSILON_F*(exp(pF)/exp(pF(+1)))^(-1-EPSILON_F)*exp(fF(+1))); //E32
exp(fF)=exp(ptilF)^(1-EPSILON_F)*exp(yF)*(EPSILON_F-1)/EPSILON_F+BETTA*THETA_F*exp(v(+1))/exp(v)*(exp(lam(+1))/exp(lam)*(exp(pic)^VARTHETA_F*pi_ss^(1-VARTHETA_F)/exp(pic(+1)))^(1-EPSILON_F)*(exp(ptilF)/exp(ptilF(+1)))^(1-EPSILON_F)*(exp(pF)/exp(pF(+1)))^-EPSILON_F*exp(fF(+1))); //E33
exp(fW)=exp(wtil)^-EPSILON_W*exp(hd)*exp(mcW)+BETTA*THETA_W*exp(v(+1))/exp(v)*(exp(lam(+1))/exp(lam)*(exp(pic)^VARTHETA_W*pi_ss^(1-VARTHETA_W)/exp(pic(+1)))^-EPSILON_W*(exp(wtil)/exp(wtil(+1)))^-EPSILON_W*(exp(w)/exp(w(+1)))^(-1-EPSILON_W)*exp(fW(+1))); //E34
exp(fW)=exp(wtil)^(1-EPSILON_W)*exp(hd)*(EPSILON_W-1)/EPSILON_W+BETTA*THETA_W*exp(v(+1))/exp(v)*(exp(lam(+1))/exp(lam)*(exp(pic)^VARTHETA_W*pi_ss^(1-VARTHETA_W)/exp(pic(+1)))^(1-EPSILON_W)*(exp(wtil)/exp(wtil(+1)))^(1-EPSILON_W)*(exp(w)/exp(w(+1)))^-EPSILON_W*exp(fW(+1))); //E35
1=(1-THETA_W)*exp(wtil)^(1-EPSILON_W)+THETA_W*(exp(w(-1))/exp(w)*exp(pic(-1))^VARTHETA_W*pi_ss^(1-VARTHETA_W)/exp(pic))^(1-EPSILON_W); //E36
exp(DelW)=(1-THETA_W)*exp(wtil)^-EPSILON_W+THETA_W*(exp(w(-1))/exp(w)*exp(pic(-1))^VARTHETA_W*pi_ss^(1-VARTHETA_W)/exp(pic))^-EPSILON_W; //E37
exp(h)=exp(hd)*exp(DelW); //E38

// Condiciones de equilibrio: Bancos (8)

exp(varrhoL)=BETTA/exp(a)*(exp(v(+1))/exp(v))*(exp(lam(+1))/exp(lam))*((1-ww_ss)*((exp(RL(+1))/exp(pic(+1)))-(exp(R)/exp(pic(+1))))+ww_ss*(exp(L(+1))/exp(L))*exp(a)*exp(varrhoL(+1))); //E39
exp(varrhoN)=BETTA/exp(a)*(exp(v(+1))/exp(v))*(exp(lam(+1))/exp(lam))*((1-ww_ss)*(exp(R)/exp(pic(+1)))+ww_ss*(exp(n(+1))/exp(n))*exp(a)*exp(varrhoN(+1))); //E40
exp(LEV)=exp(varrhoN)/(mu_ss-exp(varrhoL)); //41
exp(L)=exp(LEV)*exp(n); //E42
exp(d)=exp(L)-exp(n); //E43
exp(n)=(ww_ss/exp(a(-1)))*exp(fin)*(((exp(RL)/exp(pic))-(exp(R(-1))/exp(pic)))*exp(L(-1))+(exp(R(-1))/exp(pic))*exp(n(-1)))+iota_ss*n_ss; //E44
exp(RL)/exp(pic)= ((exp(rK) + exp(q)*(1-DELTA))/(alphaLK*exp(q(-1))) + (alphaLK-1)*exp(q)*exp(k)*exp(a(-1))/(alphaLK*exp(q(-1))*exp(k(-1)))); //E45
exp(L)=alphaLK*exp(q)*exp(k); //E46

// Definiciones (18)

stb=tb/(exp(pY)*exp(y)); //E52
sdstar=exp(rer)*dstar/(exp(pY)*exp(y)); //E53
exp(gam_Y)=exp(y)/exp(y(-1))*exp(a(-1)); //E54
exp(gam_L)=exp(L)/exp(L(-1))*exp(a(-1)); //E55
exp(gam_C)=exp(c)/exp(c(-1))*exp(a(-1)); //E56
exp(gam_I)=exp(inv)/exp(inv(-1))*exp(a(-1)); //E57
exp(gam_W)=exp(w)/exp(w(-1))*exp(a(-1)); //E58       
exp(gam_K)=exp(k)/exp(k(-1))*exp(a(-1)); //E59    
exp(gam_YC)=exp(yC)/exp(yC(-1))*exp(a(-1)); //E60  
exp(gam_YF)=exp(yF)/exp(yF(-1))*exp(a(-1)); //E61 
exp(gam_YH)=exp(yH)/exp(yH(-1))*exp(a(-1)); //E62 
exp(gam_XF)=exp(xF)/exp(xF(-1))*exp(a(-1)); //E63 
exp(gam_XH)=exp(xH)/exp(xH(-1))*exp(a(-1)); //E64 
exp(gam_XHstar)=exp(xHstar)/exp(xHstar(-1))*exp(a(-1)); //E65
exp(gam_M)=exp(m)/exp(m(-1))*exp(a(-1)); //E66
exp(gam_LAM)=exp(lam)/exp(lam(-1))/exp(a(-1)); //E67
exp(piH)=exp(pH)/exp(pH(-1))*exp(pic); //E68
exp(piF)=exp(pF)/exp(pF(-1))*exp(pic); //E69

// Procesos exógenos AR(1) (13)

v-log(v_ss)=RHO_v*(v(-1)-log(v_ss))+eps_v; //E70
u-log(u_ss)=RHO_u*(u(-1)-log(u_ss))+eps_u; //E71
z-log(z_ss)=RHO_z*(z(-1)-log(z_ss))+eps_z; //E72
a-log(a_ss)=RHO_a*(a(-1)-log(a_ss))+eps_a; //E73
yCo-log(yCo_ss)=RHO_yCo*(yCo(-1)-log(yCo_ss))+eps_yCo; //E74
g-log(g_ss)=RHO_g*(g(-1)-log(g_ss))+eps_g; //E75
zeta-log(zeta_ss)=RHO_zeta*(zeta(-1)-log(zeta_ss))+eps_zeta; //E76
eR-log(eR_ss)=RHO_eR*(eR(-1)-log(eR_ss))+eps_eR; //E77
Rstar-log(Rstar_ss)=RHO_Rstar*(Rstar(-1)-log(Rstar_ss))+eps_Rstar; //E78
ystar-log(ystar_ss)=RHO_ystar*(ystar(-1)-log(ystar_ss))+eps_ystar; //E79
pistar-log(pistar_ss)=RHO_pistar*(pistar(-1)-log(pistar_ss))+eps_pistar; //E80
pCostar-log(pCostar_ss)=RHO_pCostar*(pCostar(-1)-log(pCostar_ss))+eps_pCostar; //E81
fin-log(fin_ss)=RHO_fin*(fin(-1)-log(fin_ss))-eps_fin; //E73

// Función de Utilidad y definición de Bienestar (Consumidor y Bancos) (3)

Util=exp(v)*(log(exp(c)-VARSIGMA*exp(c(-1))) - KAPPA*((exp(h)^(1+PHI))/(1+PHI))); //E83
Welfare=Util + BETTA*Welfare(+1); //E84
Welfare_bank=exp(varrhoL)*exp(L) + exp(varrhoN)*exp(n); //E85

end;

%----------------------------------------------------------------
% 3. Estado estacionario
%----------------------------------------------------------------

steady_state_model; 

// Computando estado estacionario y parámetros endógenos
BETTA=a_ss*pi_ss/R_ss;
RL_ss=pi_ss*(GAMMA_ss+R_ss/pi_ss);
piS_ss=a_ss*pi_ss/(BETTA*Rstar_ss*xi_ss);
pistar_ss=pi_ss/piS_ss;
ptilH_ss=1;
ptilF_ss=1;
wtil_ss=1;
DelH_ss=ptilH_ss^-EPSILON_H;
DelF_ss=ptilF_ss^-EPSILON_F;
DelW_ss=wtil_ss^-EPSILON_W;
mcH_ss=(EPSILON_H-1)/EPSILON_H*ptilH_ss;
mcF_ss=(EPSILON_F-1)/EPSILON_F*ptilF_ss;
mcW_ss=(EPSILON_W-1)/EPSILON_W*wtil_ss;
hd_ss=h_ss/DelW_ss;
fW_ss=wtil_ss^-EPSILON_W*hd_ss*mcW_ss/(1-BETTA*THETA_W);
q_ss=1/u_ss;
rK_ss=q_ss*(alphaLK*(RL_ss/pi_ss)-1+DELTA-(alphaLK-1)*a_ss);
w_ss=(ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)*pH_ss*mcH_ss*z_ss*a_ss^(1-ALPHA)/rK_ss^ALPHA)^(1/(1-ALPHA));
k_ss=a_ss*ALPHA*w_ss*h_ss/((1-ALPHA)*rK_ss);
yH_ss=z_ss*(k_ss/a_ss)^ALPHA*(a_ss*h_ss)^(1-ALPHA)/DelH_ss;
fH_ss=ptilH_ss^-EPSILON_H*yH_ss*mcH_ss/(1-BETTA*THETA_H);
inv_ss=k_ss*(1-(1-DELTA)/a_ss)/u_ss;
pF_ss=(1/O-(1-O)/O*pH_ss^(1-ETA))^(1/(1-ETA));
rer_ss=mcF_ss*pF_ss;
pYy_ss=pH_ss*yH_ss/(1-sCo_ss-pF_ss*(1-mcF_ss*DelF_ss)*O*pF_ss^-ETA*(1-stb_ss));
tb_ss=stb_ss*pYy_ss;
g_ss=sg_ss*pYy_ss;
yCo_ss=sCo_ss*pYy_ss/(rer_ss*pCostar_ss);
xHstar_ss=yH_ss-(1-O)*pH_ss^-ETA*(pYy_ss-tb_ss);
xH_ss=yH_ss-xHstar_ss;
xF_ss=(pH_ss*xHstar_ss+rer_ss*pCostar_ss*yCo_ss-tb_ss)/rer_ss;
yF_ss=xF_ss;
fF_ss=ptilF_ss^-EPSILON_F*yF_ss*mcF_ss/(1-BETTA*THETA_F);
m_ss=yF_ss*DelF_ss;
yC_ss=(xF_ss/O)*pF_ss^ETA;
c_ss=yC_ss-g_ss-inv_ss;
y_ss=c_ss+inv_ss+g_ss+xHstar_ss+yCo_ss-m_ss;
pY_ss=(c_ss+inv_ss+g_ss+tb_ss)/y_ss;
lam_ss=(c_ss-VARSIGMA*c_ss/a_ss)^(-SIGMA)-BETTA*VARSIGMA*(c_ss*a_ss-VARSIGMA*c_ss)^(-SIGMA);
KAPPA=mcW_ss*lam_ss*w_ss/h_ss^PHI;
dstar_ss=-(tb_ss-(1-CHI)*rer_ss*pCostar_ss*yCo_ss)/(rer_ss*(1-Rstar_ss*xi_ss/pistar_ss/a_ss));
OSTAR=(xHstar_ss/ystar_ss)*(pH_ss/rer_ss)^ETASTAR;
sdstar_ss=rer_ss*dstar_ss/pYy_ss;
//-------------------------------------------------------------------------
ww_ss=a_ss*(1-iota_ss)/((RL_ss/pi_ss-R_ss/pi_ss)*LEV_ss+R_ss/pi_ss);
varrhoN_ss=BETTA*(1-ww_ss)*(R_ss/pi_ss)/(a_ss*(1-BETTA*ww_ss));
varrhoL_ss=BETTA*(1-ww_ss)*(RL_ss/pi_ss-R_ss/pi_ss)/(a_ss*(1-BETTA*ww_ss));
mu_ss = varrhoL_ss + varrhoN_ss/LEV_ss;
L_ss=alphaLK*q_ss*k_ss;
n_ss=L_ss/LEV_ss;
d_ss=L_ss-n_ss;
//-------------------------------------------------------------------------
Util_ss=v_ss*(log(c_ss-VARSIGMA*c_ss)-KAPPA*(h_ss^(1+PHI))/(1+PHI));
Welfare_ss=Util_ss/(1-BETTA);
Welfare_bank_ss=varrhoL_ss*L_ss + varrhoN_ss*n_ss;

// Valores iniciales para el solver numérico
lam=log(lam_ss);
c=log(c_ss);
h=log(h_ss);
hd=log(hd_ss);
w=log(w_ss);
wtil=log(wtil_ss);
mcW=log(mcW_ss);
fW=log(fW_ss);
DelW=log(DelW_ss);
inv=log(inv_ss);
k=log(k_ss);
rK=log(rK_ss);
q=log(q_ss);
y=log(y_ss);
yC=log(yC_ss);
yF=log(yF_ss);
yH=log(yH_ss);
xF=log(xF_ss);
xH=log(xH_ss);
xHstar=log(xHstar_ss);
R=log(R_ss);
xi=log(xi_ss);
RL=log(RL_ss);
pic=log(pi_ss);
rer=log(rer_ss);
pH=log(pH_ss);
ptilH=log(ptilH_ss);
pF=log(pF_ss);
ptilF=log(ptilF_ss);
pY=log(pY_ss);
piS=log(piS_ss);
mcH=log(mcH_ss);
fH=log(fH_ss);
DelH=log(DelH_ss);
mcF=log(mcF_ss);
fF=log(fF_ss);
DelF=log(DelF_ss);
dstar=   (dstar_ss);
m=log(m_ss);
tb=   (tb_ss);
v=log(v_ss);
u=log(u_ss);
z=log(z_ss);
a=log(a_ss);
zeta=log(zeta_ss);
eR=log(eR_ss);
Rstar=log(Rstar_ss);
pistar=log(pistar_ss);
pCostar=log(pCostar_ss);
yCo=log(yCo_ss);
ystar=log(ystar_ss);
g=log(g_ss);
fin=log(fin_ss);
stb=   (stb_ss);
sdstar=   (sdstar_ss);
gam_Y=log(a_ss);
gam_L=log(a_ss);
gam_C=log(a_ss);
gam_I=log(a_ss);
gam_W=log(a_ss);  
gam_K=log(a_ss);
gam_YC=log(a_ss);
gam_YF=log(a_ss);
gam_YH=log(a_ss);
gam_XF=log(a_ss);
gam_XH=log(a_ss);
gam_XHstar=log(a_ss);
gam_M=log(a_ss);
gam_LAM=log(1/a_ss);
piH=log(pi_ss);
piF=log(pi_ss);
varrhoL=log(varrhoL_ss);
varrhoN=log(varrhoN_ss);
LEV=log(LEV_ss);
L=log(L_ss);
d=log(d_ss);
n=log(n_ss);
Util=Util_ss;
Welfare=Welfare_ss;
Welfare_bank=Welfare_bank_ss;
end;

%----------------------------------------------------------------
% 4. Cómputo
%----------------------------------------------------------------

resid(1);
//options_.noprint=1;
steady;
check;

shocks;

// Choques estructurales
var eps_v; stderr 0; //SIG_v;
var eps_z; stderr 0; //SIG_z;
var eps_a; stderr 0; //SIG_a;
var eps_yCo; stderr 0; //SIG_yCo;
var eps_g; stderr 0; //SIG_g;
var eps_zeta; stderr 0; //SIG_zeta;
var eps_eR; stderr 0; //SIG_eR;
var eps_Rstar; stderr 0; //SIG_Rstar;
var eps_ystar; stderr 0; //SIG_ystar;
var eps_pistar; stderr 0; //SIG_pistar;
var eps_pCostar; stderr 0; //SIG_pCostar;
var eps_u; stderr 0; //SIG_u;
var eps_fin; stderr SIG_fin;
end;
                

stoch_simul(order=1, irf=20,noprint)y c inv rer R pic L n LEV q tb ;      