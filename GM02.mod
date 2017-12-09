%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - NKE %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo NKE basico loglinealizado. 
// El Estado Estacionario es cero debido a la loglinealizacion manual.
// Modelo presentado en Gali y Monacelli (2005).
// Politica monetaria bajo inflacion domestica.
// (c) Carlos Rojas Quiroz 

var pih x y ynat rnat r s pi p ph e ystar pistar n nx real_wage a c deprec_rate ;
varexo eps_star eps_a;
parameters sigma eta gamma phi epsilon theta beta alpha phi_pi rhoa rhoy omega
sigma_a Theta lambda kappa_a Gamma Psi;
sigma       = 1;
eta         = 1 ;
gamma       = 1;
phi         = 3;
epsilon     = 6;
theta       = 0.75;
beta        = 0.99;
alpha       = 0.4;
phi_pi      = 1.5;
rhoa        = 0.9;                                                                
rhoy        = 0.86;  
rho         = beta^(-1)-1;
omega       = sigma*gamma+(1-alpha)*(sigma*eta-1);
sigma_a     = sigma/((1-alpha)+alpha*omega);
Theta       = (sigma*gamma-1)+(1-alpha)*(sigma*eta-1);
lambda      = (1-(beta*theta))*(1-theta)/theta;
kappa_a     = lambda*(sigma_a+phi);
Gamma       = (1+phi)/(sigma_a+phi);
Psi         = -Theta*sigma_a/(sigma_a+phi);

model(linear);
x    = x(+1) - sigma_a^(-1)*(r - pih(+1) - rnat) ;                              
pih  = beta * pih(+1)+ kappa_a*x;                                                
rnat = -sigma_a*Gamma*(1-rhoa)*a + alpha*sigma_a*(Theta+Psi)*(ystar(+1)-ystar);
ynat = Gamma*a + alpha*Psi*ystar;                                                 
x    = y - ynat;                                                               
y = ystar + sigma_a^(-1)*s;
pi   = pih + alpha*(s-s(-1));
s    = s(-1) + e - e(-1) + pistar - pih;
pistar = 0;
y = a + n;
nx = alpha*(omega/sigma-1)*s;
y = c+alpha*omega/sigma*s;
real_wage = sigma*c+phi*n;
a    = rhoa*a(-1) + eps_a;
ystar= rhoy*ystar(-1) + eps_star;
r = phi_pi*pi;
pi   = p - p(-1);
pih  = ph - ph(-1);
deprec_rate=e-e(-1);
end;

shocks;
var eps_a = 1;
end;

stoch_simul(order=1,irf=20,nograph);
