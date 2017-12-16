function [RHO,sigma_lomega, mu, Eomega_down, ALPHA_omega, GAMMA_omega, CE_K, WE_K, nu, a1, b1, a3, b3, a4, b4, a12, b12, a14, b14] = BGG_ss(beta_C,RK,delta, alpha, B_K, gamma, DEF_RATE);

%Solving parameters for omega_bar in financial contract given B_K, DEF_RATE, RHO

% beta_C  = 0.99;
% RK      = ((1/beta_C)^4+0.03)^(1/4);
% delta   = 0.025;
% alpha   = 0.35;
% B_K     = 0.5;
% gamma   = 1-0.025;
RHO     = RK*beta_C;
% DEF_RATE = 0.0075;
x0(1) = log(0.3);
d_aux = 0.0001;
iter2 = 1;
tol2  = 10;
while tol2>1e-20 & iter2<200
    Fx0  = fun_bgg_ss(x0, B_K, DEF_RATE, RHO);
    dFx0 = (fun_bgg_ss(x0+d_aux,B_K, DEF_RATE, RHO)-fun_bgg_ss(x0, B_K, DEF_RATE, RHO))/d_aux;
    x1   = x0 - Fx0/dFx0;
    tol2  = abs(Fx0);
    x0   = x1;
    iter2= iter2+1;
end
omega_bar = exp(x0(1));
z = norminv(DEF_RATE);
sigma_lomega = (z^2-2*log(omega_bar))^0.5 +z;
F_omega = normcdf((log(omega_bar)+0.5*sigma_lomega^2)/sigma_lomega,0,1); %phi_omega
f_omega= exp(-0.5*(log(omega_bar)+0.5*sigma_lomega^2)^2/sigma_lomega^2)/(omega_bar*sigma_lomega*(2*pi)^0.5); %dphi_omega 
Eomega_up= normcdf((0.5*sigma_lomega^2-log(omega_bar))/sigma_lomega,0,1);
Eomega_down = 1-normcdf((0.5*sigma_lomega^2-log(omega_bar))/sigma_lomega,0,1);

mu = (RHO-1)/(Eomega_down + (omega_bar*f_omega/(1-F_omega))*(Eomega_up-omega_bar*(1-F_omega)));

ALPHA_omega= Eomega_up - omega_bar*(1-F_omega); %f_omega
GAMMA_omega= omega_bar*(1-F_omega)+ (1-mu)*Eomega_down; %g_omega
dALPHA_omega= -(1-F_omega); %df_omega
dGAMMA_omega= 1-F_omega - mu*omega_bar*f_omega; %dg_omega
rho_omega  = 1/(GAMMA_omega-ALPHA_omega*dGAMMA_omega/dALPHA_omega);
B_K_omega  = GAMMA_omega*rho_omega;
K_N        = 1/(1-B_K_omega);
df_omega   = -exp(-0.5*(log(omega_bar)+0.5*sigma_lomega^2)^2/sigma_lomega^2)/(omega_bar^2*sigma_lomega*(2*pi)^0.5)*(1+(log(omega_bar)+0.5*sigma_lomega^2)/sigma_lomega^2); %d2phi_omega
d2ALPHA_omega  = f_omega; %d2f_omega
d2GAMMA_omega  = -f_omega*(1-mu*(log(omega_bar)+0.5*sigma_lomega^2)/sigma_lomega^2); %d2g_omega
drho_omega = -1/(rho_omega^2)*(ALPHA_omega*dGAMMA_omega*d2ALPHA_omega/(dALPHA_omega^2)-ALPHA_omega*d2GAMMA_omega/dALPHA_omega);
drho_dk    = 1/(K_N^2)/(GAMMA_omega+rho_omega*dGAMMA_omega/drho_omega);
N_K        = 1-B_K_omega;
nu         = drho_dk/N_K/rho_omega;
WE_K       = N_K-gamma*ALPHA_omega*rho_omega/beta_C;
Y_K        = (RK-(1-delta))/alpha;
CE_K       = (1-gamma)*ALPHA_omega*rho_omega/beta_C;
CE_Y       = CE_K/Y_K;
dEomega_down_omega = exp(-0.5*(-log(omega_bar)+0.5*sigma_lomega^2)^2/sigma_lomega^2)/(omega_bar*sigma_lomega*(2*pi)^0.5);

delta_sigma  = 0.00001;
sigma_lomega_aux = sigma_lomega-delta_sigma;
F_omega_aux  = normcdf((log(omega_bar)+0.5*sigma_lomega_aux^2)/sigma_lomega_aux,0,1);
f_omega_aux  = exp(-0.5*(log(omega_bar)+0.5*sigma_lomega_aux^2)^2/sigma_lomega_aux^2)/(omega_bar*sigma_lomega_aux*(2*pi)^0.5);
Eomega_down_aux = 1-normcdf((0.5*sigma_lomega_aux^2-log(omega_bar))/sigma_lomega_aux,0,1);
GAMMA_omega_aux = omega_bar*(1-F_omega_aux)+ (1-mu)*Eomega_down_aux;
dGAMMA_sigma    = (GAMMA_omega-GAMMA_omega_aux)/delta_sigma;

Eomega_up_aux   = normcdf((0.5*sigma_lomega_aux^2-log(omega_bar))/sigma_lomega_aux,0,1);
ALPHA_omega_aux = Eomega_up_aux - omega_bar*(1-F_omega_aux);
dALPHA_sigma    = (ALPHA_omega-ALPHA_omega_aux)/delta_sigma;

dALPHA_omega_aux = -(1-F_omega_aux);
dGAMMA_omega_aux = 1-F_omega_aux - mu*omega_bar*f_omega_aux;
rho_omega_aux   = 1/(GAMMA_omega_aux-ALPHA_omega_aux*dGAMMA_omega_aux/dALPHA_omega_aux);
drho_sigma      = (rho_omega-rho_omega_aux)/delta_sigma;

dEomega_down_sigma = (Eomega_down-Eomega_down_aux)/delta_sigma;

a1  = dEomega_down_omega/Eomega_down*omega_bar;
b1  = dEomega_down_sigma/Eomega_down*sigma_lomega;

a3  = dALPHA_omega/ALPHA_omega*omega_bar;
b3  = dALPHA_sigma/ALPHA_omega*sigma_lomega;

a4  = drho_omega/rho_omega*omega_bar;
b4  = drho_sigma/rho_omega*sigma_lomega;

a12 = a3;
b12 = b3;

a14 = dGAMMA_omega/GAMMA_omega*omega_bar;
b14 = dGAMMA_sigma/GAMMA_omega*sigma_lomega;


