%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC01b.mod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','z'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes','Capital','Productividad'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
dates = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_z';
size_shock = 1;
resp_mat0 = [];
for ii=1:nvar
    eval(['y1=',char(varble(ii)),'_',char(shock),';']);
    y1= y1*size_shock;
    y1=[0;y1];
    resp_mat0 = [resp_mat0 y1];
end
save('Model00','resp_mat0');

load Model00

figure(1)
for ii=1:nvar
    subplot(2,4,ii);
    %plot(dates,resp_mat0(:,ii)*100,'-r','LineWidth',1.5); 
    m=plot(dates,resp_mat0(:,ii)*100);
    set(m,'Color',[51/255 130/255 214/255],'LineWidth',2.5,'LineStyle','-');
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
        legend('Modelo RBC basico');
   end
end

%%

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC02.mod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','z'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes','Capital','Productividad'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
dates = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_z';
size_shock = 1;
resp_mat1 = [];
for ii=1:nvar
    eval(['y1=',char(varble(ii)),'_',char(shock),';']);
    y1= y1*size_shock;
    y1=[0;y1];
    resp_mat1 = [resp_mat1 y1];
end
save('Model01','resp_mat1');

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC02b.mod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','z'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes','Capital','Productividad'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
dates = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_z';
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

figure(2)
for ii=1:nvar
    subplot(2,4,ii);
    plot(dates,resp_mat1(:,ii)*100,'-r',dates,resp_mat2(:,ii)*100,'--b','LineWidth',1.5); 
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
        legend('L variable','L fijo');
   end
end
