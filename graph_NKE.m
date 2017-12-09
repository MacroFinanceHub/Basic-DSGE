%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare GM01.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'pih','x','pi','s','e','r','ph','p'};
names  ={'inflacion domestica','Brecha de Producto','Inflacion CPI','Terminos de intercambio','Tipo de cambio nominal','Tasa de interes nominal', 'Precio domestico', 'IPC'};
[nper,junk1] = size(y_eps_a);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'eps_a';
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
dynare GM02.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'pih','x','pi','s','e','r','ph','p'};
names  ={'inflacion domestica','Brecha de Producto','Inflacion CPI','Terminos de intercambio','Tipo de cambio nominal','Tasa de interes nominal', 'Precio domestico', 'IPC'};
[nper,junk1] = size(y_eps_a);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'eps_a';
size_shock = 1;
resp_mat2 = [];
for ii=1:nvar
    eval(['y2=',char(varble(ii)),'_',char(shock),';']);
    y2= y2*size_shock;
    y2=[0;y2];
    resp_mat2 = [resp_mat2 y2];
end
save('Model02','resp_mat2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare GM03.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'pih','x','pi','s','e','r','ph','p'};
names  ={'inflacion domestica','Brecha de Producto','Inflacion CPI','Terminos de intercambio','Tipo de cambio nominal','Tasa de interes nominal', 'Precio domestico', 'IPC'};
[nper,junk1] = size(y_eps_a);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'eps_a';
size_shock = 1;
resp_mat3 = [];
for ii=1:nvar
    eval(['y3=',char(varble(ii)),'_',char(shock),';']);
    y3= y3*size_shock;
    y3=[0;y3];
    resp_mat3 = [resp_mat3 y3];
end
save('Model03','resp_mat3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare GM04.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'pih','x','pi','s','e','r','ph','p'};
names  ={'inflacion domestica','Brecha de Producto','Inflacion CPI','Terminos de intercambio','Tipo de cambio nominal','Tasa de interes nominal', 'Precio domestico', 'IPC'};
[nper,junk1] = size(y_eps_a);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'eps_a';
size_shock = 1;
resp_mat4 = [];
for ii=1:nvar
    eval(['y4=',char(varble(ii)),'_',char(shock),';']);
    y4= y4*size_shock;
    y4=[0;y4];
    resp_mat4 = [resp_mat4 y4];
end
save('Model04','resp_mat4');

load Model01
load Model02
load Model03
load Model04

figure(1)
for ii=1:nvar
    subplot(4,2,ii);
    plot(fechas,resp_mat1(:,ii),'--g',fechas,resp_mat2(:,ii),'-*r',fechas,resp_mat3(:,ii),'-.',fechas,resp_mat4(:,ii),'-b','LineWidth',1.5); 
    grid on; xlim([1 20]);
    hold on; 
    plot([0 20],[0 0],'-k','LineWidth',1.5)
    hold off;
    if ii>8
        xlabel('Trimestres','Fontsize',8)
    end
    ylabel('Desv. % EE','Fontsize',8)
    title(names(ii),'Interpreter','none','Fontsize',10);
   if ii==nvar
        legend('DITR','CITR','PEG','Optimal');
   end
end