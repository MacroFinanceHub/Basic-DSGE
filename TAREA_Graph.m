%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBCTAREA.mod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','lambda'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes','Capital','Choque'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_lambda';
size_shock = 1;
resp_mat0 = [];
for ii=1:nvar
    eval(['y1=',char(varble(ii)),'_',char(shock),';']);
    y1= y1*size_shock;
    y1=[0;y1];
    resp_mat0 = [resp_mat0 y1];
end
save('Model00','resp_mat0');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBCTAREA.mod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','gamma'};
names  ={'PBI','Consumo', 'Inversion','Empleo','Salario real','Tasa de interes','Capital','Choque'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
fechas = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_gamma';
size_shock = -1;
resp_mat1 = [];
for ii=1:nvar
    eval(['y1=',char(varble(ii)),'_',char(shock),';']);
    y1= y1*size_shock;
    y1=[0;y1];
    resp_mat1 = [resp_mat1 y1];
end
save('Model01','resp_mat1');

%%
load Model00
load Model01

figure(1)
for ii=1:nvar
    subplot(2,4,ii);
    hold on;
    %plot(dates,resp_mat0(:,ii)*100,'-r','LineWidth',1.5);
    %plot(dates,resp_mat1(:,ii)*100,'-b','LineWidth',1.5);
    m=plot(fechas,resp_mat0(:,ii)*100);
    set(m,'Color',[51/255 130/255 214/255],'LineWidth',2.5,'LineStyle','-');
    n=plot(fechas,resp_mat1(:,ii)*100);
    set(n,'Color',[51/255 130/255 214/255],'LineWidth',2.5,'LineStyle','--');
    hold off;
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
        legend('\lambda_t', '\gamma_t');
   end
end
