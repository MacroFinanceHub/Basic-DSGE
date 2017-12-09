%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC14.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','inom','p','pic','M'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes real','Tasa de interes nominal', 'Precios', 'Inflacion','Saldos nominales'};
[nper,junk1] = size(y_e_m);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_m';
size_shock = 1;
resp_mat1 = [];
for ii=1:nvar
    eval(['y1=',char(varble(ii)),'_',char(shock),';']);
    y1= y1*size_shock;
    y1=[0;y1];
    resp_mat1 = [resp_mat1 y1];
end
save('Model01','resp_mat1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC15.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','inom','p','pic','M'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes real','Tasa de interes nominal', 'Precios', 'Inflacion','Saldos nominales'};
[nper,junk1] = size(y_e_m);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_m';
size_shock = 1;
resp_mat2 = [];
for ii=1:nvar
    eval(['y2=',char(varble(ii)),'_',char(shock),';']);
    y2= y2*size_shock;
    y2=[0;y2];
    resp_mat2 = [resp_mat2 y2];
end
save('Model02','resp_mat2');

load Model01
load Model02

figure(1)
for ii=1:nvar
    subplot(2,5,ii);
    plot(fechas,resp_mat1(:,ii),'-r',fechas,resp_mat2(:,ii),'--b','LineWidth',1.5); 
    grid on; xlim([1 40]);
    hold on; 
    plot([0 40],[0 0],'-k','LineWidth',1.5)
    hold off;
    if ii>8
        xlabel('Trimestres','Fontsize',8)
    end
    ylabel('Desv. % EE','Fontsize',8)
    title(names(ii),'Interpreter','none','Fontsize',10);
   if ii==nvar
        legend('CIA','MIU');
   end
end