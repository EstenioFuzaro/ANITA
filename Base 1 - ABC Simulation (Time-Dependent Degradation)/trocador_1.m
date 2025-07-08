%% ------------------------------------------------------------------------- %%
% Trocador de calor concêntrico - contra-corrente                             %
%                                                                             %
%                         GEOMETRIA DA ALGETEC                                %
%                                                                             %
% Vitória Batista Godoy e Estênio Fuzaro de Almeida                           %
% -------------------------------------------------------------------------- %%
%% INICIALIZAÇÃO
clearvars; close all; clc;

nfonte = 24;
marcador=15;


mySeed = 160403;

rng_stream = RandStream('mt19937ar', 'Seed', mySeed);
RandStream.setGlobalStream(rng_stream);
%---------------------------------------------------------------------------
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');
%--------------------------------------------------------------------------
%% CALIBRAÇÃO ABC

%Trocador contra corrente
% Código de calibração bayesiana ABC -
% Considerando decaimento de U em 1 anos - U = Ui*exp(-x/a)
% Ui como variável aleatória
%a => taxa de incrustação

% Schematic representation
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tsf <-------------------------------------------< Tef : temp. casco
% ____________________________________
% Teq  >------------------------------------------->  Tsq : temp. tubo
% ____________________________________
% Tsf <-------------------------------------------< Tef : temp. casco
%
% Tsf : Temperatura de saida frio
% Tef : Temperatura de entrada frio
% Tsq : Temperatura de saida quente
% Teq : Temperatura de entrada quente
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARÂMETROS DO TROCADOR
% 
% Comprimento = 1 m
% 	Diâmetro externo do casco = 28.15 mm
% 	Diâmetro externo do tubo = 15.20 mm
% 	Espessura das paredes = 2 mm
% 	Diâmetro interno do casco (De - t) = 26.15 mm
% 	Diâmetro interno do casco (De - t) = 13.20 mm
% 
% Condições de contorno
% 	Fluxo de massa no inlet do tubo = 6.99E-2 kg/s
% 	Fluxo de massa no inlet do casco = 6.00E-2 kg/s
% 	Temperatura do fluido quente (tubo) = 90 °C
% 	Temperatura do fluido frio (casco) = 20 °C

% Trocador de calor
L       =   1.0;                % m
D_c     =   0.02815;            % m
D_t     =   0.01520;            % m
Dh      =   sqrt(D_c^2-D_t^2);  % m
Area    =   pi*D_t*L;           % m^2
Vol_f   =   pi*L*(Dh/2)^2;      % m^3
Vol_q   =   pi*L*(D_t/2)^2;     % m^3
T_ef    =   20;                 % °C
T_eq    =   90;                 % °C

%--------------------------------------------------------------------------
% Propriedades do fluido
rho_aq   =   965.3;
rho_af   =   998;     
Vzq     =   6.99E-2;        % kg/s
Vzf     =   6.00E-2;        % kg/s

k_aq     =   0.675;        % W/m.K   
cp_aq    =   4206;         % J/Kg.K
mu_aq    =   3.15E-4;      % Pa.s    

k_af     =   0.598;        % W/m.K   
cp_af    =   4182;         % J/Kg.K
mu_af    =   10.02E-4;      % Pa.s 

Pr_q    =   (mu_aq * cp_aq) / k_aq;          
Pr_f    =   (mu_af * cp_af) / k_af;        

%--------------------------------------------------------------------------
% Reynolds

Re_f    =   (4*Vzf)/(mu_af*Dh*pi);
Re_q    =   (4*Vzq)/(mu_aq*D_t*pi);
Nu_f    =   0.023*Re_f^(4/5)*Pr_f^0.4;
Nu_q    =   0.023*Re_q^(4/5)*Pr_q^0.3;
h_f     =   Nu_f*(k_af/Dh);
h_q     =   Nu_q*(k_aq/D_t);
U       =   ((1/h_q)+(1/h_f))^-1;
Umed    =   U;

%% INCRUSTAÇÃO
% 1 anos => 2190s + 400s (transiente)
% Coeficiente de incrustação => Uf=ax+b
transiente = 4000; % pontos

%--------------------------------------------------------------------------
casco = struct();
tubo = struct();
% Setting integration parameters
tubo.Te = T_eq ;  % Temperatura entrada - quente 
casco.Te = T_ef ; % Temperature entrada - saida
anos = 4;
if anos == 1
    corte = 0;
end
if anos == 4
    corte = 365*24;
end
passo = 0.1;
tt= (anos*365*24 + transiente)*passo;
intt = 0:passo:tt;           % tempo de integração 
y0 = [ tubo.Te, casco.Te ,Umed]; % condições iniciais


% components of the problem
b1=Vzq/(Vol_q*rho_aq); %1a equação - tubo 
b2=Vzf/(Vol_f*rho_af); %2a equação - casco


a1 =-Area/(Vol_q*rho_aq*cp_aq);    %1a equação tubo
a2 = Area/(Vol_f*rho_af*cp_af);    %1a equação tubo

a = 4000;   %fator de incrustação
odefun = @(t) (Umed/a) * -exp(-t / a);

opts = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[t,S]=ode45(@(t,T)funcao_contracorrente(t,T,a1,b1,a2,b2,T_eq,T_ef,odefun(t)),[intt],[y0],opts);

r1=t(transiente+corte:end); 

tubo.Ts_D = S(transiente+corte:end,1); % Temperatura saida quente
casco.Ts_D = S(transiente+corte:end,2); % Temperatura saida frio
Uf_D = S(:,3);

tubo.Ts_RP_D =tubo.Ts_D(end);   %Temperatura de saída tubo em regime permanente 
casco.Ts_RP_D =casco.Ts_D(end); %Temperatura de saída casco em regime permanente 

% 
r1 = linspace(0, anos*12, length(casco.Ts_D));
figure
plot(Uf_D)

%% Curvas de temperatura - deterministico

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.6 0.55]);
box on
hold on
set(gca, 'FontSize', 15)
%title('Fluido quente - Tubo','FontSize',13)
xlabel('Months','FontSize',24)
ylabel('Temperature ($^\circ$C)','FontSize',24)  %Arrumar latex
xlim([0,anos*12])
plot(r1,tubo.Ts_D,'-r','LineWidth',1.7)

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.6 0.55]);
box on
hold on
set(gca, 'FontSize', 15)
xlabel('Months','FontSize',24)
ylabel('Temperature ($^\circ$C)','FontSize',24)
%yticks(0:1:24)
xlim([0,anos*12])
plot(r1,casco.Ts_D,'-b','LineWidth',1.7)
%legend('','FontSize',20)

figure 
plot(t,S(:,3))

%% Estocástico 

% Generating training data
n = 200;     % número de dados experimentais gerados
g = 0;
k = 625;        % Shape parameter for Gamma distribution
theta = Umed / k;  % Scale parameter for Gamma distribution

for g = 1:n
    Ue(g) = gamrnd(k, theta);  % U inicial sem incrustação
    a1 = -Area / (Vol_q * rho_aq * cp_aq);    % 1a equação tubo
    a2 = Area / (Vol_f * rho_af * cp_af);     % 2a equação casco
    y0 = [tubo.Te, casco.Te, Ue(g)];  % condições iniciais
    a = 4000;  % fator de incrustação
    odefun = @(t) (Ue(g) / a) * -exp(-t / a);

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
    [t, S] = ode45(@(t, T) funcao_contracorrente(t, T, a1, b1, a2, b2, T_eq, T_ef, odefun(t)), [intt], [y0], opts);

    r1 = t(transiente + corte:end);

    tubo.Ts_e(g, :) = S(transiente + corte:end, 1);  % Curva Temperatura fluido quente - tubo
    casco.Ts_e(g, :) = S(transiente + corte:end, 2);  % Curva Temperatura fluido frio - casco
    U_e(g, :) = S(transiente + corte:end, 3);  % Curva do coef. de transf. de calor U

    tubo.Ts_RP_e(g) = tubo.Ts_e(g, end);  % Temperatura de saída em regime permanente
    casco.Ts_RP_e(g) = casco.Ts_e(g, end);  % Temperatura de saída em regime permanente
    Uf_e(g) = U_e(g, end);
end
g = 0;

%% Dados para validação 
% Generating validation data
n = 10;     % número de dados experimentais gerados
g = 0;

for g = 1:n
    Uv(g) = gamrnd(k, theta);  % U inicial sem incrustação
    a = 4000;  % fator de incrustação
    odefun = @(t) (Uv(g) / a) * -exp(-t / a);
    a1 = -Area / (Vol_q * rho_aq * cp_aq);  % 1a equação tubo
    a2 = Area / (Vol_f * rho_af * cp_af);   % 2a equação casco
    y0 = [tubo.Te, casco.Te, Uv(g)];  % condições iniciais

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
    [t, S] = ode45(@(t, T) funcao_contracorrente(t, T, a1, b1, a2, b2, T_eq, T_ef, odefun(t)), [intt], [y0], opts);

    r1 = t(transiente + corte:end);

    tubo.Ts_v(g, :) = S(transiente + corte:end, 1);  % Curva Temperatura fluido quente - tubo
    casco.Ts_v(g, :) = S(transiente + corte:end, 2);  % Curva Temperatura fluido frio - casco
    U_v(g, :) = S(transiente + corte:end, 3);  % Curva do coef. de transf. de calor U

    tubo.Ts_RP_v(g) = tubo.Ts_v(g, end);  % Temperatura de saída em 4 anos
    casco.Ts_RP_v(g) = casco.Ts_v(g, end);  % Temperatura de saída em 4 anos
end
g = 0;



%% Plot Distribuição experimental
Nbins = round(sqrt(n)); 
[x_bins,x_freq,x_area] = randvar_pdf(Uf_e,Nbins);
%                       = randvar(dados, Nbins)
[x_ksd_e,x_supp1_e] = ksdensity(Uf_e);

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.6 0.5]);
% [left bottom width height]
box on
hold on
h(1) = bar(x_bins,x_freq,1.0);
h(2) = plot(x_supp1_e,x_ksd_e,'-k','LineWidth',1.5);
%xline(a, '--r', 'LineWidth', 3); 
hold on; grid on;
set(gca, 'FontSize', 15)
legend('Experimental',"Mean",...
    'FontSize', 16);
xlim([min(Uf_e)*0.9,max(Uf_e)*1.1])
xlabel('$\mathbf{U}$ [$\frac{W}{m^{2}K}$]','FontSize',16);
ylabel('Frequency','FontSize',16);


%% Dados para validação 
% Generating validation data
n = 10;     % número de dados experimentais gerados
g = 0;
k = 625;        % Shape parameter for Gamma distribution
theta = Umed / k;  % Scale parameter for Gamma distribution

for g = 1:n
    Uv(g) = gamrnd(k, theta);  % U inicial sem incrustação
    a = 4000;  % fator de incrustação
    odefun = @(t) (Uv(g) / a) * -exp(-t / a);
    a1 = -Area / (Vol_q * rho_aq * cp_aq);  % 1a equação tubo
    a2 = Area / (Vol_f * rho_af * cp_af);   % 2a equação casco
    y0 = [tubo.Te, casco.Te, Uv(g)];  % condições iniciais

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
    [t, S] = ode45(@(t, T) funcao_contracorrente(t, T, a1, b1, a2, b2, T_eq, T_ef, odefun(t)), [intt], [y0], opts);

    r1 = t(transiente + corte:end);

    tubo.Ts_v(g, :) = S(transiente + corte:end, 1);  % Curva Temperatura fluido quente - tubo
    casco.Ts_v(g, :) = S(transiente + corte:end, 2);  % Curva Temperatura fluido frio - casco
    U_v(g, :) = S(transiente + corte:end, 3);  % Curva do coef. de transf. de calor U

    tubo.Ts_RP_v(g) = tubo.Ts_v(g, end);  % Temperatura de saída em 4 anos
    casco.Ts_RP_v(g) = casco.Ts_v(g, end);  % Temperatura de saída em 4 anos
end
g = 0;

 
% %% Calibração 
% % 2- PRIOR DISTRIBUTION
% %__________________________________________________________________________
% % % Set bound vectors:
Umin = Umed*0.85; % minimum values;
Umax = Umed*1.15; % maximum values.



%% Inicialização 
%Numero de simulações
NMC=1000;
%Tolerância 
tol=0.002;
%Erro inicial (alto)
erro=1e6;

%ABC Parameters 
Parameters.Ts_tubo = tubo.Ts_RP_D;
Parameters.Texp_tubo = tubo.Ts_RP_e;
Parameters.Ts_casco = casco.Ts_RP_D;
Parameters.Texp_casco = casco.Ts_RP_e;

% For now, this is the minimum error value:
JM=fun_objetiva(Parameters)
XMLE = U;
YX_casco = casco.Ts_D;
YX_tubo = tubo.Ts_D;
U=250; % Chute inicial do U incrustado
% At first, the number of accepted samples is null:
Naccept = 0;
tmp_erro = [];
%% Calibração
%% 
% MONTE CARLO SIMULATION:
tic
for ii = 1:NMC
        U_prop= unifrnd(Umin,Umax);      % U inicial sem incrustação 
        Utent(ii) = U_prop;
        %a = (-U_prop/2)*(1/2500);
        odefun = @(t) (U_prop/a) * -exp(-t / a);
        a1=-Area/(Vol_q*rho_aq*cp_aq);    %1a equação tubo
        a2= Area/(Vol_f*rho_af*cp_af);    %1a equação tubo
        y0 = [ tubo.Te, casco.Te ,U_prop]; % condições iniciais


        opts = odeset('RelTol',1e-6,'AbsTol',1e-6); 
        [t,S]=ode45(@(t,T)funcao_contracorrente(t,T,a1,b1,a2,b2,T_eq,T_ef,odefun(t)),[intt],[y0],opts);

        r1=t(transiente+corte:end);

        tubo.Ts = S(transiente+corte:end,1);               %Curva Temperatura fluido quente - tubo
        casco.Ts = S(transiente+corte:end,2);              %Curva Temperatura fluido frio - casco
        Ut(ii,:) = S(transiente+corte:end,3);              %Curva do parâmetro U 

        tubo.Ts_RP = tubo.Ts(end);     %Temperatura de saída em 4 anos
        casco.Ts_RP = casco.Ts(end);   %Temperatura de saída em 4 anos



        Parameters.Ts_tubo = tubo.Ts_RP;     %Temperatura de saída em 4 anos
        Parameters.Ts_casco = casco.Ts_RP;   %Temperatura de saída em 4 anos
        Uf = Ut(ii,end);                     %Valores de U em 4 anos
        Uf_tent(ii) = Uf;

        [JT] = fun_objetiva(Parameters); 
        err(ii)=JT;
        % X,YX and PX get candidate values if the it is accepted:
            if JT<=tol
            U = Uf;
            YX_casco = casco.Ts;
            YX_tubo = tubo.Ts;
            J = JT;
            Naccept = Naccept+1;
            tmp_erro = [JT; tmp_erro];

            % For minimum error :
            % if(J < tol)
            %     XMLE = U;
            % end
        end

        % Markov chain gets X, YX and PX. In the the case the current candidate
        % is not accepted, these values are already for the last accepted
        % candidate:
        samples(ii,:) = U;             %Candidatos aceitos
        Y_casco(ii,:) = YX_casco;      %Resposta dos candidatos aceitos
        Y_tubo(ii,:) = YX_tubo;        %Resposta dos candidatos aceitos
        % Jsamples(ii,1) = J;            %Erro

        %Print its values for monitoring:
        %fprintf('Acceptance rate: %f\n', Naccept/ii); % 40-50 % --> Adjust random walk step size.
        %fprintf('%.5f completed\n', ii/NMC);
end
toc


%% 5.1.1 - TRACE PLOTS
parstr = {'$k$'};
Umin=min(samples)
Umax=max(samples)

%% 5.1.2 - DEFINE BURN-IN AND GET SAMPLES

burnin = 0.0*NMC;            % 20% burnin samples;
x = samples(burnin+1:end,:); % Get random variables samples;
U = U(burnin+1:end,1);       % Get the generated displacements
Ns = size(x,1);              % Get number of random variariable semples.

%% 5.1.3 - PLOT SAMPLES

for par = 1:length(Umin)
    xmean = mean(x(:,par)); % Compute sample mean;
    xstd = std(x(:,par));   % Compute sample standard deviation;
    Pc = 95;                % Define 95th percentile
    r_p = 0.5*(100+Pc); x_upp = prctile(x(:,par),r_p); % Get upper;
    r_m = 0.5*(100-Pc); x_low = prctile(x(:,par),r_m); % Get lower percentile.

    figure
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.5 0.55]);
    % [left bottom width height]
    h(1) = plot(x(:,par),'oc');
    hold on
    h(2) = plot([1 Ns],[xmean xmean],'--k','linewidth',1.4);
    h(3) = plot([1 Ns],[xmean-xstd xmean-xstd],'-.b','linewidth',1.4);
    h(4) = plot([1 Ns],[x_low x_low],':m','linewidth',1.4);
    h(5) = plot([1 Ns],[mean(Uf_e) mean(Uf_e)],'-r','linewidth',1.4);
    plot([1 Ns],[xmean+xstd xmean+xstd],'-.b','linewidth',1.4)
    plot([1 Ns],[x_upp x_upp],':m','linewidth',1.4)
    grid on   
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    xlabel('Samples','Interpreter','latex','fontsize',16)
    ylabel(parstr{par},'Interpreter','latex','fontsize',16)
    legstr = {'Samples','Mean','Mean $\pm$ std','$95\%$ interval', 'Real Value'};
    legend(h([1,2,3,4,5]),legstr,'Interpreter','latex','fontsize',15)
end

%% 5.1.4 - CONVERGENCE STUDY: 1ST MOMENT ESTIMATOR PLOT

for par = 1:length(Umin)
    xcmean = cumsum(x(:,par))'./(1:Ns);
    figure
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.5 0.55]);
% [left bottom width height]
    h(1) = plot(x(:,par));
    hold on
    h(2) = plot(xcmean,'--k','linewidth',1.0); hold on
    h(3) = plot(Umin(par)*ones(1,Ns),'--r','linewidth',1.0);
    plot(Umax(par)*ones(1,Ns),'--r','linewidth',1.0)
    grid on 
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    xlabel('N','Interpreter','latex','fontsize',1)
    ylabel('$\mathbf{U}$ [$\frac{W}{m^{2}K}$]','Interpreter','latex','fontsize',16)
    legstr = {'Samples','Mean','Support'};
    legend(h([1,2,3]),legstr,'Interpreter','latex')
    %title('1st central moment','Interpreter','latex')
end
%% 5.1.6 - PLOT PARAMETERS DISTRIBUTIONS
% 
 Nbins = round(sqrt(Ns));
 Nksd = round(sqrt(Ns));
 for par = 1:length(Umin)

     a = min(samples)*0.9; b = max(samples)*1.1;
     [x_bins,x_freq,x_area] = randvar_pdf(x(:,par),Nbins);
     [x_ksd,x_supp1] = ksdensity(x(:,par),'Support',[a b]);
     [x_cdf,x_supp2] = ecdf(x(:,par));

     figure
     set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.6 0.5]);
% [left bottom width height]
    h(1) = bar(x_bins,x_freq,1.0);
    hold on
    h(2) = plot(x_supp1,x_ksd,'-k','LineWidth',1.5);
    h(3) = line([a b],[1/(b-a) 1/(b-a)],'Color','r','LineStyle','--','LineWidth',1.0);
    h(4) = plot(x_supp1_e,x_ksd_e,'-g','LineWidth',1.5);                              %Experimental
    line([a a],[0 1/(b-a)],'Color','r','LineStyle','--','LineWidth',1.5)
    line([b b],[0 1/(b-a)],'Color','r','LineStyle','--','LineWidth',1.5)
    xlim([Umin(par)-(Umax(par)-Umin(par))*(0.2),Umax(par)+(Umax(par)-Umin(par))*(0.2)])
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    xlabel('$\mathbf{U}$ [$\frac{W}{m^{2}K}$]','Interpreter','latex','fontsize',16)
    ylabel('Density','Interpreter','latex','fontsize',16)%Mudar
    legend(h([1,2,3,4]),'Chain 1','EPDF','Prior','Target','Interpreter','latex')
 end
    %% 

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.6 0.5]);
% [left bottom width height]
plot(Uf_tent,'+b','MarkerSize',6)
hold on
plot(samples,'.r','MarkerSize',8) 
grid on
xlabel('N','Interpreter','latex','fontsize',1)
ylabel('$\mathbf{U}$ [$\frac{W}{m^{2}K}$]','Interpreter','latex')
%legend(h([3,1,2]),'Prior','Chain 1','EPDF','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',15)

%% Priori
Nbins= round(sqrt(NMC));
[x_bins,x_freq,x_area] = randvar_pdf(Uf_tent,Nbins);
%                       = randvar(dados, Nbins)

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.5 0.5]);
% [left bottom width height]
box on
hold on
h(1) = bar(x_bins,x_freq,1.0);
% Armazena o objeto de histograma em 'h'
hold on; grid on;
set(gca, 'FontSize', 15)
%legend('Priori',...
%    'FontSize', 20, 'Location', 'NorthWest');
xlabel('$\mathbf{U}$ [$\frac{W}{m^{2}K}$]','FontSize',16);
%xlim([Umin(par)-(Umax(par)-Umin(par))*(0.1),Umax(par)+(Umax(par)-Umin(par))*(0.1)])
ylabel('Frequency','FontSize',16);

%% PLots

%Definindo pontos aleatórios de validação para serem plotados
num_colunas_plotar = 50;         % Quantidade de colunas aleatórias para plotar por linha
num_colunas_total = size(tubo.Ts_v, 2);
indices_colunas = zeros(size(tubo.Ts_v, 1), num_colunas_plotar);
for i = 1:size(tubo.Ts_v, 1)
    indices_colunas(i, :) = randperm(num_colunas_total, num_colunas_plotar);
end

mq=mean(Y_tubo);
mf=mean(Y_casco);

r1 = linspace(0, anos*12, length(casco.Ts_D));
r1=r1';

%limites do envelope de probabilidade
limsup_q=prctile(Y_tubo,99);
liminf_q=prctile(Y_tubo,1);

limsup_f=prctile(Y_casco,99);
liminf_f=prctile(Y_casco,1);

x_confq =[r1', fliplr(r1')];
y_confq=[limsup_q, fliplr(liminf_q)];

x_conff =[r1', fliplr(r1')];
y_conff=[limsup_f, fliplr(liminf_f)];


%Curvas de temperatura com intervalo de confiança
figure
ax1 = subplot(2,1,1);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.05 0.5 0.8]);
% [left bottom width height]
box on
hold on
set(gca, 'FontSize', 15)
title('Hot fluid - Tube','FontSize',16)
xlabel('Time[Months]','FontSize',16)
ylabel('Temperature [$^\circ$C]','FontSize',16)  %Arrumar latex
yticks(0:1:90)
fill(x_confq, y_confq, 1,'Facecolor', "#FF0000",'EdgeColor','none','facealpha', 0.3); hold on
hold on
plot(r1,mq,'--r','LineWidth',1.5)
hold on
for i = 1:size(tubo.Ts_v, 1)  % Para cada linha em A
    plot(indices_colunas(i, :)*(anos*12/length(r1)), tubo.Ts_v(i, indices_colunas(i, :)), 'ok','MarkerSize',5);  % 'o-' para marcar os pontos com círculos e conectar com linhas
end
hold on
xlim([0,anos*12])
legend('99$\%$ Interval',"Mean",'Validation','FontSize',15,'Location','southeast')
%saveas(gcf,'intervalo_q1','epsc')

%figure
ax1 = subplot(2,1,2);
box on
hold on
set(gca, 'FontSize', 15)
title('Cold fluid - Shell','FontSize',16)
xlabel('Time[Months]','FontSize',16)
ylabel('Temperatura [$^\circ$C]','FontSize',16)
fill(x_conff, y_conff, 1,'Facecolor', "#3141f5",'EdgeColor','none','facealpha', 0.3); hold on
hold on
plot(r1,mf,'--b','LineWidth',1.5) 
hold on
for i = 1:size(tubo.Ts_v, 1)  % Para cada linha em A
    plot(indices_colunas(i, :)*(anos*12/length(r1)), casco.Ts_v(i, indices_colunas(i, :)), 'ok','MarkerSize',5);  % 'o-' para marcar os pontos com círculos e conectar com linhas
end
hold on
xlim([0,anos*12])
legend('$99\%$ Interval',"Mean",'Validation','FontSize',15)
%saveas(gcf,'intervalo_f1','epsc')

%% VARIÁVEIS PARA TRANSFER LEARNING - HE1
month_tube_he1 = [];
temp_tube_he1 = [];
month_shell_he1 = [];
temp_shell_he1 = [];
time_vector_he1 = []; % Novo vetor de tempo para as amostras

for i = 1:size(tubo.Ts_v, 1)
    % Gerar os meses e selecionar as temperaturas amostradas
    month_tube = indices_colunas(i, :)*(anos*12/length(r1));
    temp_tube = tubo.Ts_v(i, indices_colunas(i, :));
    month_shell = indices_colunas(i, :)*(anos*12/length(r1));
    temp_shell = casco.Ts_v(i, indices_colunas(i, :));
    
    % Concatenar nos arrays
    month_tube_he1 = [month_tube_he1 month_tube];
    temp_tube_he1 = [temp_tube_he1 temp_tube];
    month_shell_he1 = [month_shell_he1 month_shell];
    temp_shell_he1 = [temp_shell_he1 temp_shell];
    
    % Salvar o vetor de tempo amostrado
    time_vector_he1 = [time_vector_he1 r1(indices_colunas(i, :))];
end

% Garantir que o vetor de tempo seja linearizado (caso não esteja)
time_vector_he1 = time_vector_he1(:)';

% Remover valores duplicados em time_vector_he1
[time_vector_he1_unique, unique_indices] = unique(time_vector_he1);

% Atualizar as temperaturas para corresponder aos índices únicos
noisy_temp_tube_he1_unique = temp_tube_he1(unique_indices);
noisy_temp_shell_he1_unique = temp_shell_he1(unique_indices);

% Adicionar ruído às temperaturas (apenas após a remoção de duplicatas)
std_temp_tube_he1 = std(noisy_temp_tube_he1_unique);
noise_tube = randn(size(noisy_temp_tube_he1_unique));
scaled_noise_tube = 0.1 * std_temp_tube_he1 * noise_tube;
noisy_temp_tube_he1 = noisy_temp_tube_he1_unique + scaled_noise_tube;

std_temp_shell_he1 = std(noisy_temp_shell_he1_unique);
noise_shell = randn(size(noisy_temp_shell_he1_unique));
scaled_noise_shell = 0.1 * std_temp_shell_he1 * noise_shell;
noisy_temp_shell_he1 = noisy_temp_shell_he1_unique + scaled_noise_shell;

% Interpolação para Vetor de Tempo Regularizado
% Obter o intervalo mínimo e máximo de tempo no vetor único
time_min = min(time_vector_he1_unique); % Tempo inicial (mínimo global)
time_max = max(time_vector_he1_unique); % Tempo final (máximo global)

% Criar o vetor de tempo igualmente espaçado (500 pontos)
num_points = 500; % Sabemos que são 500 pontos no total
time_vector_regular = linspace(time_min, time_max, num_points);

% Interpolação das temperaturas para o vetor de tempo regularizado
noisy_temp_tube_he1_regular = interp1(time_vector_he1_unique, noisy_temp_tube_he1, time_vector_regular, 'linear');
noisy_temp_shell_he1_regular = interp1(time_vector_he1_unique, noisy_temp_shell_he1, time_vector_regular, 'linear');

% Salvar dados para 1 ano ou 4 anos, incluindo o vetor de tempo regularizado
if anos == 1
    noisy_temp_shell_he1_1year = noisy_temp_shell_he1_regular;
    noisy_temp_tube_he1_1year = noisy_temp_tube_he1_regular;
    time_vector_1year = time_vector_regular;
    save he1_data_1year noisy_temp_shell_he1_1year noisy_temp_tube_he1_1year time_vector_1year
end
if anos == 4
    noisy_temp_shell_he1_4year = noisy_temp_shell_he1_regular;
    noisy_temp_tube_he1_4year = noisy_temp_tube_he1_regular;
    time_vector_4year = time_vector_regular;
    save he1_data_4year noisy_temp_shell_he1_4year noisy_temp_tube_he1_4year time_vector_4year
end

%%

figure
set(gcf,'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.5 0.5]);
% [left bottom width height]
box on
plot(Y_casco(:,1),Y_tubo(:,1),'*')
hold on
plot(Y_casco(:,(length(r1)/2)),Y_tubo(:,(length(r1)/2)),'*g')
hold on
plot(Y_casco(:,length(r1)),Y_tubo(:,length(r1)),'*r')
xlabel('Temperature - Cold fluid [$^\circ$C]','FontSize',16)
ylabel('Temperature - Hot fluid [$^\circ$C]','FontSize',16)

% yticks(0:1:90)
% xticks(0:1:90)
legend('M\^es 0',"M\^es 24",'M\^es 48','FontSize',15)

%% Gráfico Cascata

n_tempos = anos*12;
resol = floor(length(r1)/n_tempos);
for i = 1:n_tempos
    j = i*resol - (resol-1);
    [x_ksd_e,x_supp1_e] = ksdensity(Y_tubo(:,j));
    x_values(i,:) = x_supp1_e ;
    pdf_temperatura(:,i) = x_ksd_e;
    tempos(i) = r1(j);
    temperaturas(i) = mq(j);
end
% Plote as curvas em cascata
figure
for i = 1:n_tempos
   plot3(tempos(i) * ones(1, length(x_values)), x_values(i,:), pdf_temperatura(:, i), 'r');
   hold on;
end
xlabel('Time [months]');
ylabel('Temperature [$^\circ$C]');
zlabel('Temperature PDF');
%title('PDF de Temperatura do Fluido Quente em Função do Tempo (Cascata)');
grid on;