function q = fun_bgg_ss(x, B_K, DEF_RATE, RHO)
%
%  

omega_bar = exp(x(1));

z = norminv(DEF_RATE);
sigma_lomega = (z^2-2*log(omega_bar))^0.5 +z;

phi_omega = normcdf((log(omega_bar)+0.5*sigma_lomega^2)/sigma_lomega,0,1);
dphi_omega= exp(-0.5*(log(omega_bar)+0.5*sigma_lomega^2)^2/sigma_lomega^2)/(omega_bar*sigma_lomega*(2*pi)^0.5);
Eomega_up= normcdf((0.5*sigma_lomega^2-log(omega_bar))/sigma_lomega,0,1);
Eomega_down = 1-normcdf((0.5*sigma_lomega^2-log(omega_bar))/sigma_lomega,0,1);

mu = (RHO-1)/(Eomega_down + (omega_bar*dphi_omega/(1-phi_omega))*(Eomega_up-omega_bar*(1-phi_omega)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q(1) = real((B_K-(omega_bar*(1-phi_omega)+(1-mu)*Eomega_down)*RHO)^2*100);


end

