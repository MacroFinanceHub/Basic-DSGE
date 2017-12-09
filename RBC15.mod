%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmica.
// Se incluye dinero en la funcion de utilidad.
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv z g inom p pic M m;
predetermined_variables kap;
varexo e_z e_g e_m;
parameters alpha delta betta theta psi eta rho_z rho_g rho_m 
z_ss lab_ss r_ss  kap_ss w_ss y_ss c_ss inv_ss g_ss i_ss p_ss pic_ss 
m_ss M_ss C_Y I_Y G_Y;

alpha  = 1-0.33;
delta  = 0.023;
betta  = 0.99;
theta  = 1/2.75;
eta    = 2.00;
rho_z  = 0.95;
rho_g  = 0.75;
rho_m  = 0.00;
z_ss   = 1;
G_Y    = 0.155;
i_ss   = 1/betta-1;
lab_ss = 1/((1-theta)/(alpha*theta*z_ss)*((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)+1);
y_ss   = z_ss*(((1-alpha)*betta/(1-betta+betta*delta))^((1-alpha)/alpha))*lab_ss;
w_ss   = alpha*y_ss/lab_ss;
kap_ss = (1-alpha)*betta/(1-betta+betta*delta)*y_ss;
inv_ss = delta*kap_ss;
r_ss   = (1-alpha)*y_ss/kap_ss-delta;
c_ss   = ((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)*y_ss;
M_ss   = 0.728149533266271;
p_ss   = 1;
m_ss   = log(M_ss)-log(p_ss);
psi    = (M_ss/p_ss)^(1/eta)/(c_ss*(1+i_ss)/i_ss);
pic_ss = 0;
g_ss   = G_Y*y_ss;
C_Y    = c_ss/y_ss;
I_Y    = inv_ss/y_ss;

model;
theta/exp(c)                = (1-theta)/((1-exp(lab))*exp(w));
1/exp(c)                    = betta*1/exp(c(+1))*(1+exp(r(+1)));
1+exp(r(+1))                = (1+exp(inom))/(1+pic(+1));
exp(M)/exp(p)               = psi^eta*exp(c)^eta*((1+exp(inom))/exp(inom))^eta;
m                           = M-p;
pic                         = p-p(-1);
exp(w)                      = alpha*exp(y)/exp(lab);
exp(r)+delta                = (1-alpha)*exp(y)/exp(kap);
exp(y)                      = exp(c)+exp(innv)+exp(g);
exp(kap(+1))                = (1-delta)*exp(kap)+exp(innv);
exp(y)                      = exp(z)*exp(kap)^(1-alpha)*exp(lab)^alpha;
z                           = (1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g                           = (1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
m-m(-1)+pic                 = (1-rho_m)*(pic_ss)  + rho_m*(m(-1)-m(-2)+pic(-1)) + e_m;
end;

steady_state_model;
lab =log(lab_ss);
c   =log(c_ss); 
w   =log(w_ss); 
r   =log(r_ss); 
y   =log(y_ss); 
kap =log(kap_ss); 
innv=log(inv_ss);
inom=log(i_ss);
p   =log(p_ss);
pic =pic_ss;
M   =log(M_ss);
m   =m_ss;
z   =log(z_ss);
g   =log(g_ss);
end;

shocks;

var e_z; stderr 0.01;
var e_g; stderr 0.01;
var e_m; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
