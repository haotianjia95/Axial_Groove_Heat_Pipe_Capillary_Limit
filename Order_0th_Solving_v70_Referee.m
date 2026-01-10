% THIS PROGRAM DETERMINS LEADING ORDER PROBLEM OF FLOW TRHOUGH A HEAT PIPE 
% WRITTEN BY HAOTIAN JIA, Copyright TRANSPORT PHENOMENA LAB, TUFTS UNIV.
%% Specify Fluid Proporties and Solver Settings
clear all; close all; clc; format compact; format  shortE

% Solver Settings ---------------------------------------------------------
N_An = 1000; % number of terms, Andrew's code
meMax = 0.005 ; % max mesh element size
meMin = 0.0025; % min mesh element size

% Fluid Proporties --------------------------------------------------------
    T = 300;
 NIST = readmatrix('NISTSaturationProperties_Water.csv');
p_sat = interp1(NIST(:,1),NIST(:,02),T)*1E+6
rho_v = interp1(NIST(:,1),NIST(:,15),T)
rho_l = interp1(NIST(:,1),NIST(:,03),T)
 mu_v = interp1(NIST(:,1),NIST(:,24),T)
 mu_l = interp1(NIST(:,1),NIST(:,12),T)
 h_lv = interp1(NIST(:,1),NIST(:,18)-NIST(:,06),T)*1E+3
sigma = interp1(NIST(:,1),NIST(:,14),T)
  c_l = interp1(NIST(:,1),NIST(:,10),T)
  c_v = interp1(NIST(:,1),NIST(:,22),T)

 %   Tmu = 200
 % mu_v = interp1(NIST(:,1),NIST(:,24),Tmu)
 % mu_l = interp1(NIST(:,1),NIST(:,12),Tmu)

%  Trho = 200
% rho_v = interp1(NIST(:,1),NIST(:,15),Trho)
% rho_l = interp1(NIST(:,1),NIST(:,03),Trho)

   ta = deg2rad(33) % Contract angle
   % ta = deg2rad(157.4 - 0.55* T)

%% Define Cross-Sectional Domain Geometry
Ng = 26; % number of grooves
tp = 2*pi/Ng;
tg = 7.9*pi/180; 

% Rg_star = Rg_starArray(rgcntr)
% Rg_star = [0.00118 0.00120 0.00122] %
Rg_star = [0.00122] %
N_geom = length(Rg_star);

RL_star = 0.001091;
Rp_star_Ctra = RL_star*sin(tg/2)/sin(pi/2 - ta + tg/2)
Rp_star_Geom = RL_star*sin(tg/2)
Rp_star = max(Rp_star_Ctra,Rp_star_Geom);
 P_star = tp*RL_star;

clearence = Rg_star - RL_star*sin(tg/2) - Rp_star 

% % Non-dimensionalization ------------------------------------------------
rho = rho_v/rho_l % Density   Ratio 
 mu = mu_v/mu_l   % Viscosity Ratio
 RL = RL_star/P_star           
 Rp = Rp_star/P_star
 RgArr = Rg_star/P_star 

% % theta_m^0 and R_m^0 array ---------------------------------------------
N_tm0 = 200;
% tm0_min = 1.3930;
% tm0_deno = sin(tm0_max+tg/2)
tm0_min =      - tg/2 + 0.001;
tm0_max = pi/2 - tg/2 - 0.100;
tm0_deno = linspace(sin(tm0_min+tg/2),sin(tm0_max+tg/2),N_tm0); % this will keep the Δp evenly spaced
tm0 = asin(tm0_deno) - tg/2;
% tm0 = [tm0(tm0<0) 0 tm0(tm0>0)];
% clear tm0_deno; tm0_deno = sin(tm0 + tg/2);

tm0_zero = logspace(-7,-2,50);  % This section used for quantify when approach zero
tm0 = [tm0(tm0<0) 0 tm0_zero tm0(tm0>1E-2)];

clear tm0_deno; tm0_deno = sin(tm0 + tg/2);

Rm0 = RL*sin(tg/2)./tm0_deno; %  greater than RL*sin(tg/2), no upper limit
xm = RL*cos(tg/2) - sqrt(Rm0.^2 - RL^2 * sin(tg/2)^2); % Rm center x

% figure
% plot(tm0,Rm0,'.')
% xline(pi/2 - tg/2,'-','pi/2 - tg/2')
% yline(Rp,':','Rp,min')
% xlabel('\theta_m^0')
% ylabel('R_m^0')

% % Creat Data Structures -------------------------------------------------
timestr = datetime('now','Format','yyyy-MM-dd');
% filename = sprintf('XiRgs%1.2i-%1.2i_T%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,timestr)
% filename = sprintf('XiRgs%1.2i-%1.2i_T%i_Tmu%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,Tmu,timestr)
filename = sprintf('XiRgs%1.2i-%1.2i_T%i_Trho%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,Trho,timestr)

%% Calculate \tilde\xi^00 & \tilde\xi^01 Values
tXi0N = NaN(N_geom,length(tm0));
tQv0N = NaN(N_geom,length(tm0));
tXi0A = NaN(N_geom,2);
tQv0A = NaN(N_geom,4);


for i_geom = 1:length(RgArr) % loop counter for different domain shape

    Rg = RgArr(i_geom)
    % Analytical Calculation ----------------------------------------------
    % initial guess analytical:
    if i_geom == 1; tXi00_ig = -180 %-48.93 + 0.8744*T - 0.005149*T^2; % -180 for 200C
    else          ; tXi00_ig = 0.8*tXi0A(i_geom-1,1);    end % 0.9 for 200C
    % initial guess analytical:
    [tXi00,tXi01] = XiRootAndrew(rho,mu,tp,tg,RL,Rg,N_An,tXi00_ig);
    [tQv00,tQl00,tQv01,tQl01] = A_FullMassFlowRate(mu,tp,tg,RL,Rg,tXi00,tXi01,N_An);
    % Numerical Calculation -----------------------------------------------
    for i = 1:length(Rm0) % loop counter for different meniscus raidus Rm
        fprintf('Geom %i/%i, Iteration %i/%i, Rm0 = %1.2i.\n',i_geom,length(RgArr),i,length(Rm0),Rm0(i));

        Model = ModelingFunction(tp,tg,RL,Rg,Rm0(i),xm(i)); % domain only determined by [tp tg RL Rg]
        if tm0(i) <= 1.3930
            Model.Mesh = generateMesh(Model,'Hmax',meMax); % meshing from DGM
        else
            Model.Mesh = generateMesh(Model,'Hmax',meMin); % meshing from DGM
        end
        
        % initial guess numerical:
        if i <= 3
            tXi0N_ig = tXi00 + tm0(i).*tXi01  ;
        else
            Np = (i < 5)*i + (i > 5)*(5);
            p0 = polyfit(tm0(1:i-1),tXi0N(i_geom,1:i-1),Np);
            tXi0N_ig = polyval(p0,tm0(i));
        end

        tXi0N(i_geom,i) = fzero(@(tXi) XiNumericalSolver(Model,rho,mu,tXi,'rootFinding'),tXi0N_ig);
        tQv0N(i_geom,i) = XiNumericalSolver(Model,rho,mu,tXi0N(i_geom,i),'tQvCalculate');

    end % end of i loop
    tXi0A(i_geom,:) = [tXi00,tXi01];
    tQv0A(i_geom,:) = [tQv00,tQl00,tQv01,tQl01];
         
end

save(filename) % Saving Xi Numerical Results

% %% Plotting Section (~\xi^0 Problem VERISION 1) -------------------------------------
% clear all; close all; clc; format  compact; format long
% 
% load('XiRgs1.22e+00-1.22e+00_T60_2024-12-28.mat'); cstr = '#1368AA';
% % load('XiRgs1.22e+00-1.22e+00_T200_2024-12-28.mat'); cstr = '#D95319';
% % load('XiRgs1.22e+00-1.22e+00_T300_2024-12-29.mat'); cstr = '#7E2F8E';
% 
% % load('XiRgs1.18e+00-1.22e+00_T200_2025-01-09.mat'); cstr = '#D95319';
% 
% ResultNumrical = tXi0N;
% ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
% ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
% 
% for ig = 1:length(tXi0A(:,1))
%     hfig = figure;
% 
%     subplot(2,1,1);
%     plot(tm0,ResultNumrical(ig,:),'.',DisplayName='$\tilde{\xi}^{0}$');hold on;
%     yline(ResultAnalytic(ig,:),'-',LineWidth=1.5,Color='red',DisplayName='$\tilde{\xi}^{0,0}$')
%     xlabel('$\theta_m^0~\mathrm{[rad]}$');
%     ylabel('a) $\tilde{\xi}$ and $\tilde{\xi}^{0}$')
%     xlim([0 inf])
%     legend(Location='southwest')
% 
%     subplot(2,1,2);
%     nonIdx = find(tm0>1E-3);
%     slope1 = ResultAbsDiff(ig,nonIdx(1))/tm0(nonIdx(1));
%     slope2 = ResultAbsDiff(ig,nonIdx(1))/tm0(nonIdx(1))^2;
%     loglog(tm0,ResultAbsDiff(ig,:),'.',Color='red',DisplayName='$|\tilde{\xi} - \tilde{\xi}^{0}|$'); hold on;
%     loglog(tm0,slope1*tm0*1.25,'-',Color='black',DisplayName='slope of 1'); hold on;
%     % loglog(tm0,slope2*tm0.^2,'--' ,Color='black',DisplayName='slope of 2'); hold on;
%     ylabel('b) $|\tilde{\xi} - \tilde{\xi}^{0}|$')
%     xlabel('$\theta_m^0~\mathrm{[rad]}$');
%     xlim([1E-6 1.6])
%     ylim([min(ResultAbsDiff(ig,:)),inf])
%     legend(Location='northwest')
% 
%     grid on
%     hfigname = ['Figures\xi^0_vs_theta_m^0' sprintf('_T%i',T)];
%     picturewidth = 17; % set this parameter and keep it forever
%     hw_ratio = 0.65; % feel free to play with this ratio
%     set(findall(hfig,'-property','FontSize'  ),'FontSize'  ,17)
%     set(findall(hfig,'-property','MarkerSize'),'MarkerSize',11)
%     % set(findall(hfig,'-property','Box'),'Box','off') % optional
%     set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
%     set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
%     set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
%     pos = get(hfig,'Position');
%     set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%     print(hfig,hfigname,'-dpdf','-painters','-fillpage')
%     print(hfig,hfigname,'-dpng','-painters')
% end

%% Plotting Section (~\xi^0 Problem VERISION 2.2, Paper Plot) -------------------------------------
% Rg_star = 1.22 mm groove depth
clear all; close all; clc; format  compact; format long

hfig = figure('DefaultAxesFontSize',11);
fsLgd = 10
fsAxi = 14

tiledlayout(4,1,"Padding","tight")

nexttile
    load('XiRgs1.22e+00-1.22e+00_T60_2024-12-28.mat');
    load('XiCaseABRgs1.18e+00-1.22e+00_T60_2025-07-17.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    ResultAnaly_CA = tXi00_CaseA(1,3);
    ResultAnaly_CB = tXi00_CaseB(1,3);
    % load('XiCaseANumRgs1.22e+00-1.22e+00_T60_2025-03-24');
    % ResultNumcl_CA = tXi00_CaseANum
    load('XiCaseBNumRgs1.22e+00-1.22e+00_T60_2025-07-17');
    ResultNumcl_CB = tXi00_CaseBNum
    caseAoffpercent =  abs((ResultAnaly_CA - tXi0N(15))/tXi0N(15))*100
    % subplot(5,1,1);
    p1 = plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='A',Color='#77AC30');hold on;
    % plot([0],ResultNumcl_CA,'s',LineWidth=1.5,DisplayName='$\tilde{\xi}^{\rm{A,Num}}$',Color='#77AC30')
    p2 = plot([0],ResultAnaly_CB,'+',LineWidth=1.5,DisplayName='B',Color='black')
    % plot([0],ResultNumcl_CB,'o',LineWidth=1.5,DisplayName='B numerical',Color='#EDB120')
    p3 = plot(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='C',Color='#0072BD');
    % plot([0],ResultAnalytic,'*',LineWidth=1.5,DisplayName='$\tilde{\xi}^{0}$')
    % xline(pi/2,'--','$\pi/2$',Interpreter='latex',FontSize=fsLgd,HandleVisibility='off',...
        % LabelVerticalAlignment='middle',LabelHorizontalAlignment='left')
    uistack(p1, 'top');
    uistack(p2, 'top');
    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-10 100]); ylim([-30 -5]); % ylim([-10^(3) 0]); %
    xticks([-10 0 10 20 30 40 50 60 70 80 90])
    xticklabels({'-10°' '0°' '10°' '20°' '30°' '40°' '50°' '60°' '70°' '80°' '90°'})
    lgd = legend(Location='northeast',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=60^{\circ}$C',FontSize=fsLgd)
  
nexttile
    load('XiRgs1.22e+00-1.22e+00_T200_2024-12-28.mat');
    load('XiCaseABRgs1.18e+00-1.22e+00_T200_2025-07-17.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    ResultAnaly_CA = tXi00_CaseA(1,3);
    ResultAnaly_CB = tXi00_CaseB(1,3);
    % load('XiCaseANumRgs1.22e+00-1.22e+00_T200_2025-03-24');
    % ResultNumcl_CA = tXi00_CaseANum
    load('XiCaseBNumRgs1.22e+00-1.22e+00_T200_2025-07-17');
    ResultNumcl_CB = tXi00_CaseBNum
    caseAoffpercent =  abs((ResultAnaly_CA - tXi0N(15))/tXi0N(15))*100
    % subplot(5,1,2);
    p1 = plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='A',Color='#77AC30');hold on;
    % plot([0],ResultNumcl_CA,'s',LineWidth=1.5,DisplayName='$\tilde{\xi}^{\rm{A,Num}}$',Color='#77AC30')
    p2 = plot([0],ResultAnaly_CB,'+',LineWidth=1.5,DisplayName='B',Color='black')
    % plot([0],ResultNumcl_CB,'o',LineWidth=1.5,DisplayName='B numerical',Color='#EDB120')
    p3 = plot(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='C',Color='#0072BD');
    % plot([0],ResultAnalytic,'*',LineWidth=1.5,DisplayName='$\tilde{\xi}^{0}$')
    % xline(pi/2,'--','$\pi/2$',Interpreter='latex',FontSize=fsLgd,HandleVisibility='off',...
        % LabelVerticalAlignment='middle',LabelHorizontalAlignment='left')
    uistack(p1, 'top');
    uistack(p2, 'top');
    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-10 100]);   ylim([-200 -50]);
    xticks([-10 0 10 20 30 40 50 60 70 80 90])
    xticklabels({'-10°' '0°' '10°' '20°' '30°' '40°' '50°' '60°' '70°' '80°' '90°'})
    lgd = legend(Location='northeast',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=200^{\circ}$C',FontSize=fsLgd)

    nexttile
    load('XiRgs1.22e+00-1.22e+00_T300_2024-12-29.mat');
    load('XiCaseABRgs1.18e+00-1.22e+00_T300_2025-07-17.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    ResultAnaly_CA = tXi00_CaseA(1,3);
    ResultAnaly_CB = tXi00_CaseB(1,3);
    % load('XiCaseANumRgs1.22e+00-1.22e+00_T300_2025-03-24');
    % ResultNumcl_CA = tXi00_CaseANum
    load('XiCaseBNumRgs1.22e+00-1.22e+00_T300_2025-07-17');
    ResultNumcl_CB = tXi00_CaseBNum
    caseAoffpercent =  abs((ResultAnaly_CA - tXi0N(15))/tXi0N(15))*100
    % subplot(5,1,3);
    p1 = plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='A',Color='#77AC30');hold on;
    % plot([0],ResultNumcl_CA,'s',LineWidth=1.5,DisplayName='$\tilde{\xi}^{\rm{A,Num}}$',Color='#77AC30')
    p2 = plot([0],ResultAnaly_CB,'+',LineWidth=1.5,DisplayName='B',Color='black')
    % plot([0],ResultNumcl_CB,'o',LineWidth=1.5,DisplayName='B numerical',Color='#EDB120')
    p3 = plot(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='C',Color='#0072BD');
    % plot([0],ResultAnalytic,'*',LineWidth=1.5,DisplayName='$\tilde{\xi}^{0}$')
    % xline(pi/2,'--','$\pi/2$',Interpreter='latex',FontSize=fsLgd,HandleVisibility='off',...
        % LabelVerticalAlignment='middle',LabelHorizontalAlignment='left')
    uistack(p1, 'top');
    uistack(p2, 'top');
    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-10 100]);   ylim([-800 0])
    xticks([-10 0 10 20 30 40 50 60 70 80 90])
    xticklabels({'-10°' '0°' '10°' '20°' '30°' '40°' '50°' '60°' '70°' '80°' '90°'})
    lgd = legend(Location='northeast',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=300^{\circ}$C',FontSize=fsLgd)

nexttile    % subplot(5,1,5);
load('XiRgs1.22e+00-1.22e+00_T60_2024-12-28.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    semilogy(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='$\tilde{\xi}~[\tilde\rho(60^{\circ}\rm{C}),\tilde\mu(60^{\circ}\rm{C})]$',Color='#D95319');hold on;
    load('XiRgs1.22e+00-1.22e+00_T60_Tmu200_2025-03-15');
    ResultNumrical = tXi0N;
    semilogy(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='$\tilde{\xi}~[\tilde\rho(60^{\circ}\rm{C}),\tilde\mu(200^{\circ}\rm{C})]$',Color='#0072BD');
    % For  fictitious fluid 
    % semilogy([0],ResultAnalytic(1,:),'*',DisplayName='$\tilde{\xi}^0~[\tilde\rho(60^{\circ}\rm{C}),\tilde\mu(200^{\circ}\rm{C})]$',Color='#77AC30')
    % For  200C fluid 
    load('XiRgs1.22e+00-1.22e+00_T60_Trho200_2025-03-15');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    semilogy(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='$\tilde{\xi}~[\tilde\rho(200^{\circ}\rm{C}),\tilde\mu(60^{\circ}\rm{C})]$',Color='#EDB120')
    % semilogy([0],ResultAnalytic(1,:),'*',DisplayName='$\tilde{\xi}^0~[\tilde\rho(200^{\circ}\rm{C}),\tilde\mu(60^{\circ}\rm{C})]$',Color='#4DBEEE')
    load('XiRgs1.22e+00-1.22e+00_T200_2024-12-28.mat');
    ResultNumrical = tXi0N;
    semilogy(rad2deg(tm0),ResultNumrical(1,:),'.',DisplayName='$\tilde{\xi}~[\tilde\rho(200^{\circ}\rm{C}),\tilde\mu(200^{\circ}\rm{C})]$',Color='#77AC30')
    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-10 120]); ylim([-10^(3) -1]); %
    xticks([-10 0 10 20 30 40 50 60 70 80 90])
    xticklabels({'-10°' '0°' '10°' '20°' '30°' '40°' '50°' '60°' '70°' '80°' '90°'})
    lgd = legend(Location='east',Box='off',Orientation='vertical',FontSize=fsLgd);
    % title(lgd,'$T_{\rm{o}}^*=60^{\circ}$C',FontSize=fsLgd)

    % common x-axix label
    xlabel('$\theta_{\mathrm{m}}$',FontSize=fsAxi);

% grid on
hfigname = ['Figures\xi_vs_theta_m_v3'];
% hfigname = ['Figures\xi_vs_theta_m_v4_zoom'];
picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 1; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'  ),'FontSize'  ,16)
set(findall(hfig,'-property','MarkerSize'),'MarkerSize',10)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[2 2 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
print(hfig,hfigname,'-dpng','-painters')

%%
hfig = figure('DefaultAxesFontSize',11);
fsLgd = 10
fsAxi = 14

    load('XiRgs1.22e+00-1.22e+00_T60_2024-12-28.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    subplot(3,1,1);
    nonIdx = find(tm0>1E-3);
    slope1 = ResultAbsDiff(1,nonIdx(1))/tm0(nonIdx(1));
    loglog(tm0,ResultAbsDiff,'.',Color='#0072BD',DisplayName='$|{\tilde \xi} - {\tilde \xi}^{0}_{\rm series}|$'); hold on;
    loglog(tm0,slope1*tm0*1.25,'-',Color='black',DisplayName='slope of 1'); hold on;
    % loglog(tm0,slope2*tm0.^2,'--' ,Color='black',DisplayName='slope of 2'); hold on;
    grid on
    ylabel('$|{\tilde \xi} - {\tilde \xi}^{0}_{\rm series}|$',FontSize=fsAxi)
    xlim([1E-6 1.6])
    ylim([min(ResultAbsDiff),inf])
    lgd = legend(Location='northwest',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=60^{\circ}$C',FontSize=fsLgd)
    
    load('XiRgs1.22e+00-1.22e+00_T200_2024-12-28.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    subplot(3,1,2);
    nonIdx = find(tm0>1E-3);
    slope1 = ResultAbsDiff(1,nonIdx(1))/tm0(nonIdx(1));
    loglog(tm0,ResultAbsDiff,'.',Color='#0072BD',DisplayName='$|{\tilde \xi} - {\tilde \xi}^{0}_{\rm series}|$'); hold on;
    loglog(tm0,slope1*tm0*1.25,'-',Color='black',DisplayName='slope of 1'); hold on;
    % loglog(tm0,slope2*tm0.^2,'--' ,Color='black',DisplayName='slope of 2'); hold on;
    grid on
    ylabel('$|{\tilde \xi} - {\tilde \xi}^{0}_{\rm series}|$',FontSize=fsAxi)
    xlim([1E-6 1.6])
    ylim([min(ResultAbsDiff),inf])
    lgd = legend(Location='northwest',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=200^{\circ}$C',FontSize=fsLgd)
    
    load('XiRgs1.22e+00-1.22e+00_T300_2024-12-29.mat');
    ResultNumrical = tXi0N;
    ResultAnalytic = tXi0A(:,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    subplot(3,1,3);
    nonIdx = find(tm0>1E-3);
    slope1 = ResultAbsDiff(1,nonIdx(1))/tm0(nonIdx(1));
    loglog(tm0,ResultAbsDiff,'.',Color='#0072BD',DisplayName='$|{\tilde \xi} - {\tilde \xi}^{0}_{\rm series}|$'); hold on;
    loglog(tm0,slope1*tm0*1.25,'-',Color='black',DisplayName='slope of 1'); hold on;
    % loglog(tm0,slope2*tm0.^2,'--' ,Color='black',DisplayName='slope of 2'); hold on;
    grid on
    ylabel('$|{\tilde \xi} - {\tilde \xi}^{0}_{\rm series}|$',FontSize=fsAxi)
    xlim([1E-6 1.6])
    ylim([min(ResultAbsDiff),inf])
    lgd = legend(Location='northwest',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=300^{\circ}$C',FontSize=fsLgd)

    xlabel('$\theta_{\mathrm{m}}~\mathrm{[rad]}$',FontSize=fsAxi);


hfigname = ['Figures\xi_vs_theta_m_abs_diff_v2'];
picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','MarkerSize'),'MarkerSize',10)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
print(hfig,hfigname,'-dpng','-painters')


%% Plotting Section (~\xi^0 Problem VERISION 3.1) -------------------------------------
% Rg_star = all groove depth
clear all; close all; clc; format  compact; format long

ig = 1
switch ig
    case 1; hfigname = ['Figures\xi_vs_theta_m_v3_Rg118'];
    case 2; hfigname = ['Figures\xi_vs_theta_m_v3_Rg120'];
    case 3; hfigname = ['Figures\xi_vs_theta_m_v3_Rg122'];
end 
hfig = figure('DefaultAxesFontSize',11);
fsLgd = 10
fsAxi = 14

tcl = tiledlayout(3,1,"Padding","tight")

nexttile
    load('XiRgs1.18e+00-1.22e+00_T60_2025-01-07.mat');
    load('XiCaseABRgs1.18e+00-1.22e+00_T60_2025-07-17.mat');
    ResultNumrical = tXi0N(ig,:);
    ResultAnalytic = tXi0A(ig,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    ResultAnaly_CA = tXi00_CaseA(1,ig);
    ResultAnaly_CB = tXi00_CaseB(1,ig);

    % caseoffpercent =  abs((ResultAnalytic - tXi0N(ig,15))/tXi0N(ig,15))*100
    ACratioXi =  abs(ResultAnaly_CA/ResultAnalytic)
    % subplot(5,1,1);
    plot(tm0,ResultNumrical,'.',DisplayName='C',Color='#0072BD');hold on;
    plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='A',Color='#77AC30');
    plot([0],ResultAnaly_CB,'+',LineWidth=1.5,DisplayName='B',Color='#EDB120')
    
    % plot([0],ResultAnalytic,'*',LineWidth=1.5,DisplayName='C_0',Color='r')
    % plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='$\tilde{\xi}^{\rm{A}}$',Color='#EDB120');hold on;
    % xline(pi/2,'--','$\pi/2$',Interpreter='latex',FontSize=fsLgd,HandleVisibility='off',...
        % LabelVerticalAlignment='middle',LabelHorizontalAlignment='left')

    % load('XiCaseBNumRgs1.18e+00-1.18e+00_T60_2025-07-17');
    % load('XiCaseBNumRgs1.22e+00-1.22e+00_T60_2025-03-24');
    switch ig
    case 1; load('XiCaseBNumRgs1.18e+00-1.18e+00_T60_2025-07-17');
    case 2; load('XiCaseBNumRgs1.20e+00-1.20e+00_T60_2025-07-17');
    case 3; load('XiCaseBNumRgs1.22e+00-1.22e+00_T60_2025-07-17');
    end
    ResultNumcl_CB = tXi00_CaseBNum
    plot([0],ResultNumcl_CB,'o',LineWidth=1.5,DisplayName='B numerical',Color='#EDB120')

    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-0.2 1.8]); % ylim([-20 -5]); % ylim([-10^(3) 0]); %
    xticks([-0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 pi/2])
    xticklabels({'-0.2' '0' '0.2' '0.4' '0.6' '0.8' '1' '1.2' '1.4' '\pi/2'})
    lgd = legend(Location='northeast',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=60^{\circ}$C',FontSize=fsLgd)
    txt = ['$\tilde{\xi}^{\rm{A}} / \tilde{\xi}^{0}~=~$' sprintf('%1.4f',ACratioXi)];
    text(min(tm0),min(ResultNumrical)*0.7,txt,'FontSize',12,Interpreter='latex')
  
nexttile
    load('XiRgs1.18e+00-1.22e+00_T200_2025-01-09.mat');
    load('XiCaseABRgs1.18e+00-1.22e+00_T200_2025-07-17.mat');
    ResultNumrical = tXi0N(ig,:);
    ResultAnalytic = tXi0A(ig,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    ResultAnaly_CA = tXi00_CaseA(1,ig);
    ResultAnaly_CB = tXi00_CaseB(1,ig);
    % caseAoffpercent =  abs((ResultAnalytic - tXi0N(ig,15))/tXi0N(ig,15))*100
    ACratioXi =  abs(ResultAnaly_CA/ResultAnalytic)
    % subplot(5,1,1);
    plot(tm0,ResultNumrical,'.',DisplayName='C',Color='#0072BD');hold on;
    plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='A',Color='#77AC30');
    plot([0],ResultAnaly_CB,'+',LineWidth=1.5,DisplayName='B',Color='#EDB120');
    % plot([0],ResultAnalytic,'*',LineWidth=1.5,DisplayName='$\tilde{\xi}^{0}$',Color='r')
    % xline(pi/2,'--','$\pi/2$',Interpreter='latex',FontSize=fsLgd,HandleVisibility='off',...
        % LabelVerticalAlignment='middle',LabelHorizontalAlignment='left')

    switch ig
    case 1; load('XiCaseBNumRgs1.18e+00-1.18e+00_T200_2025-07-17');
    case 2; load('XiCaseBNumRgs1.20e+00-1.20e+00_T200_2025-07-17');
    case 3; load('XiCaseBNumRgs1.22e+00-1.22e+00_T200_2025-07-17');
    end
    ResultNumcl_CB = tXi00_CaseBNum
    plot([0],ResultNumcl_CB,'o',LineWidth=1.5,DisplayName='B numerical',Color='#EDB120')

    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-0.2 1.8]); % ylim([-20 -5]); % ylim([-10^(3) 0]); %
    xticks([-0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 pi/2])
    xticklabels({'-0.2' '0' '0.2' '0.4' '0.6' '0.8' '1' '1.2' '1.4' '\pi/2'})
    lgd = legend(Location='northeast',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=200^{\circ}$C',FontSize=fsLgd)
    txt = ['$\tilde{\xi}^{\rm{A}} / \tilde{\xi}^{0}~=~$' sprintf('%1.4f',ACratioXi)];
    text(min(tm0),min(ResultNumrical)*0.7,txt,'FontSize',12,Interpreter='latex')

nexttile 
    load('XiRgs1.18e+00-1.22e+00_T300_2025-01-07.mat');
    load('XiCaseABRgs1.18e+00-1.22e+00_T300_2025-07-17.mat');
    ResultNumrical = tXi0N(ig,:);
    ResultAnalytic = tXi0A(ig,1); %+ tXi0A(:,2).*tm0;
    ResultAbsDiff  = abs(ResultNumrical - ResultAnalytic);
    ResultAnaly_CA = tXi00_CaseA(1,ig);
    ResultAnaly_CB = tXi00_CaseB(1,ig);
    % caseAoffpercent =  abs((ResultAnalytic - tXi0N(ig,15))/tXi0N(ig,15))*100
    ACratioXi =  abs(ResultAnaly_CA/ResultAnalytic)
    % subplot(5,1,1);
    plot(tm0,ResultNumrical,'.',DisplayName='C',Color='#0072BD');hold on;
    plot([0],ResultAnaly_CA,'x',LineWidth=1.5,DisplayName='A',Color='#77AC30');
    plot([0],ResultAnaly_CB,'+',LineWidth=1.5,DisplayName='B',Color='#EDB120');
    % plot([0],ResultAnalytic,'*',LineWidth=1.5,DisplayName='$\tilde{\xi}^{0}$',Color='r')
    % xline(pi/2,'--','$\pi/2$',Interpreter='latex',FontSize=fsLgd,HandleVisibility='off',...
        % LabelVerticalAlignment='middle',LabelHorizontalAlignment='left')

    switch ig
    case 1; load('XiCaseBNumRgs1.18e+00-1.18e+00_T300_2025-07-17');
    case 2; load('XiCaseBNumRgs1.20e+00-1.20e+00_T300_2025-07-17');
    case 3; load('XiCaseBNumRgs1.22e+00-1.22e+00_T300_2025-07-17');
    end
    ResultNumcl_CB = tXi00_CaseBNum
    plot([0],ResultNumcl_CB,'o',LineWidth=1.5,DisplayName='B numerical',Color='#EDB120')

    ylabel('$\tilde{\xi}$',Rotation=0,FontSize=fsAxi)
    xlim([-0.2 1.8]); % ylim([-20 -5]); % ylim([-10^(3) 0]); %
    xticks([-0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 pi/2])
    xticklabels({'-0.2' '0' '0.2' '0.4' '0.6' '0.8' '1' '1.2' '1.4' '\pi/2'})
    lgd = legend(Location='northeast',Box='off',Orientation='vertical',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=300^{\circ}$C',FontSize=fsLgd)
    txt = ['$\tilde{\xi}^{\rm{A}} / \tilde{\xi}^{0}~=~$' sprintf('%1.4f',ACratioXi)];
    text(min(tm0),min(ResultNumrical)*0.7,txt,'FontSize',12,Interpreter='latex')

    % common x-axix label
    xlabel('$\theta_{\mathrm{m}}~\mathrm{[rad]}$',FontSize=fsAxi);

    tclText =  sprintf('$R_{\\rm{g}}^* = %1.2f\\rm{mm}$',Rg_star(ig)*1000)
    title(tcl,tclText,Interpreter='latex')

% grid on
% hfigname = ['Figures\xi_vs_theta_m_v2'];

picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 1; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'  ),'FontSize'  ,16)
set(findall(hfig,'-property','MarkerSize'),'MarkerSize',10)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[2 2 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
% print(hfig,hfigname,'-dpng','-painters')


%% Define Longitudinal Problem Parameters& % Plotting Section (Capillary Limit VERISION 1)
clear all; close all; clc; format  compact; format long
% load('XiRgs1.18e+00-1.22e+00_T60_2025-01-07.mat')
load('XiRgs1.18e+00-1.22e+00_T200_2025-01-09.mat')
% load('XiRgs1.18e+00-1.22e+00_T300_2025-01-07.mat')

clear Model meMax meMin clearance filename timestr NIST
clear tm0_max tm0_min Rg epsilon V_sv_star V_sl_star
clear tQv00 tQv01 tQl00 tQl01 tXi00 tXi01 tXi00_ig tXi0N_ig

% Solver Settings ---------------------------------------------------------
N_pt = 200;  % numbers of interpulations points

% Cross-Sectional Variables -----------------------------------------------
Dh = 2*tg*(Rg_star.^2 - RL_star^2)./...
    (2*(Rg_star - RL_star) + tg.*(Rg_star+RL_star) );

phi_l = Ng*tg/(2*pi);

% Wick permeability per Sparrow et al. (1963)
SparrowRatio = [ 0.7   0.8   0.85  0.9   0.95]; 
SparrowValue = [18.08 16.42 15.32 14.35 15.06;...
                15.56 14.42 14.26 14.98 17.62;...
                14.52 14.33 14.91 16.36 19.17;...
                14.26 14.84 15.82 17.51 20.14;...
                14.71 16.15 17.39 19.07 21.26;...
                15.52 17.28 18.52 20.05 21.88;...
                16.33 18.16 19.34 20.71 22.28;...
                17.04 18.86 19.95 21.18 22.56];
SparrowAngel = deg2rad([5 10 15 20 30 40 50 60]);
Sparrow_fRe = interp2(SparrowRatio,SparrowAngel,SparrowValue,RL_star./Rg_star,tg)
clear SparrowRatio SparrowValue SparrowAngel

K = Dh.^2 * phi_l./(2*Sparrow_fRe)
Aw_star = pi*(Rg_star.^2 - RL_star^2);

    % \tilde \xi for Case A & B -------------------------------------------
    tXi00_CaseA = - pi*rho*RL_star^4./(8*K.*Aw_star*mu);
    
    for i_geom = 1:length(RgArr) 
        [K1_star,K2_star,K4_star,Ca_star] = A_FullMassFlowRate_CaseAB(mu,tp,tg,RL_star,RL,Rg_star(i_geom),RgArr(i_geom),P_star,tXi0A(i_geom,1),tXi0A(i_geom,2),N_An);
        tXi00_CaseB(i_geom) = (((pi*rho*RL_star^4)/(16 *Ng*mu))-K4_star*RL_star/2) / (K1_star+K2_star-Ca_star*K4_star);
        clear K1_star K2_star K4_star Ca_star
    end
    
    timestr = datetime('2025-07-17','Format','yyyy-MM-dd');
    filename = sprintf('XiCaseABRgs%1.2i-%1.2i_T%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,timestr)
    save(filename,'tXi00_CaseA','tXi00_CaseB') % Saving Xi Numerical Results
    % ---------------------------------------------------------------------

% Longitudinal Variables --------------------------------------------------
g = 9.81;
phi = 0;

% La_star = 0.100 * 5; % This has changed to 5x larger
% Le_star = 0.100;
% Lc_star = 0.095;
% La_star = linspace(0.1,1,100);
Lt_star = 0.1 + 0.1 + 0.1; % Total pipe length
% La_star = linspace(Lt_star*0.975,Lt_star*0.025,100)
La_star = 0.2975:-0.0025:0.025;

Le_star = (Lt_star-La_star)./2 ;
Lc_star = (Lt_star-La_star)./2 ;

Leff_star = Le_star/2 + La_star + Lc_star/2;

epsilon =  P_star./La_star     ; % small parameter
Dp_cmax_star = sigma/Rp_star;
   Ca_hat = Dp_cmax_star*P_star/sigma;
V_sv_star = Dp_cmax_star*P_star*epsilon/mu_v;
V_sl_star = Dp_cmax_star*P_star*epsilon/mu_l;

    hfigx = figure;
    plot(La_star,epsilon,'.')
    xlabel('La* (m)')
    ylabel('$\epsilon$')
    box on; grid on
    picturewidth = 17; % set this parameter and keep it forever
    hw_ratio = 0.65; % feel free to play with this ratio
    set(findall(hfigx,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document
        % set(findall(hfig,'-property','Box'),'Box','on') % optional
    set(findall(hfigx,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfigx,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfigx,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfigx,'Position');
    set(hfigx,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])


% Non-dimensionalization --------------------------------------------------
  Aw = Aw_star.*K/(P_star^4);
  Le = Le_star./La_star     ;
  Lc = Lc_star./La_star     ;
Bo_v = rho_v*g.*La_star*sin(phi)/Dp_cmax_star;
Bo_l = rho_l*g.*La_star*sin(phi)/Dp_cmax_star;
d_Bo = (rho_l-rho_v)*g.*La_star*sin(phi)/Dp_cmax_star; % ΔBo

Bo_prep_v = rho_v*g*cos(phi)*P_star/Dp_cmax_star
Bo_prep_l = rho_l*g*cos(phi)*P_star/Dp_cmax_star

epsl2 = epsilon.^2;


% Capillary limit calculation and plots
clc
Qv0A_C0 = zeros(N_geom,length(La_star));
Qv0N_C0 = zeros(N_geom,length(La_star));
Qv0N_C1 = zeros(N_geom,length(La_star));

qc_CaseA = zeros(N_geom,length(La_star));
qc_CaseB = zeros(N_geom,length(La_star));

for ig = 1:N_geom % loop counter for different domain shape

    SC.rho = rho;
    SC.mu = mu;
    SC.Ng = Ng;
    SC.RL = RL;
    SC.tg = tg;
    SC.Aw = Aw(ig);
    SC.Ca_hat = Ca_hat;

    tXi00 = tXi0A(ig,1);
    tXi01 = tXi0A(ig,2);
    tQv00 = tQv0A(ig,1);
    tQv01 = tQv0A(ig,3);

    for j = 1:length(La_star)  % loop counter for different pipe length La

        % Longitudinal Varaible Structure:
        SL.Le = Le(j);
        SL.Lc = Lc(j);
        SL.d_Bo = d_Bo(j);

        % Qv^0 Initial Guess ----------------------------------------------
        if j == 1
            switch T
                case  60; tQv0_igA =  1; %0.005;%0.0009;%[0.01 1]; % 60C-1
                          tQv0_igA_57 = 0.1;
                case 200; tQv0_igA =  0.1; %0.005;%0.0009;%[0.01 1]; % 200C-0.1 
                          tQv0_igA_57 = 0.01;
                case 300; tQv0_igA =  0.006; %0.005;%0.0009;%[0.01 1]; % 300C-0.01 
                          tQv0_igA_57 = 0.003;
            end 
            tQv0_igN =  tQv0_igA ;
        else
            tQv0_igA = Qv0A_C0(ig,j-1); % tilde xi 0 initial guess 1
            tQv0_igN = Qv0N_C0(ig,j-1);
        end

        Qv0A_C0(ig,j) = fzero(@(Qv0Arg) QvRootAndrew(SC,SL,tXi00,tXi01,tQv00,tQv01,Qv0Arg,'C0'),tQv0_igA);
        Qv0N_C0(ig,j) = fzero(@(Qv0Arg) ...
            QvRootNumerical(SC,SL,tm0,Rm0,tXi0N(ig,:),tQv0N(ig,:),Qv0Arg,N_pt,'C0'),tQv0_igN);
        Qv0N_C0_57(ig,j) = fzero(@(Qv0Arg) ...
            QvRootNumerical(SC,SL,tm0,Rm0,tXi0N(ig,:),tQv0N(ig,:),Qv0Arg,N_pt,'C0_57'),tQv0_igA_57);
        Qv0N_C1(ig,j) = fzero(@(Qv0Arg) ...
            QvRootNumerical(SC,SL,tm0,Rm0,tXi0N(ig,:),tQv0N(ig,:),Qv0Arg,N_pt,'C1'),tQv0_igN);
    end

    clear SC SL

    % Case A & B
    
    qc_CaseA(ig,:) = (Dp_cmax_star - (rho_v - rho_l)*g*sin(phi)*Lt_star)./...
    (8*mu_v*Leff_star/(pi*rho_v*RL_star^4*h_lv)+...
    mu_l*Leff_star./(K(ig)*rho_l*Aw_star(ig)*h_lv));

    [K1_star,K2_star,K4_star,Ca_star] = A_FullMassFlowRate_CaseAB(mu,tp,tg,RL_star,RL,Rg_star(ig),RgArr(ig),P_star,tXi0A(ig,1),tXi0A(ig,2),N_An);
    qc_CaseB(ig,:) = (Dp_cmax_star - (rho_v - rho_l)*g*sin(phi)*Lt_star)./...
    (Leff_star*( 8*mu_v/(pi*rho_v*RL_star^4*h_lv) + (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv) )/(K1_star + K2_star - Ca_star*K4_star)*mu_l ));
    % qc_CaseB(ig,:) = (Dp_cmax_star - (rho_v - rho_l)*g*sin(phi)*Lt_star)./...
    % (Leff_star*( 8*mu_v/(pi*rho_v*RL_star^4*h_lv) + (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv) )/(K1_star + K2_star)*mu_l ));


end
tQv0N(:,1);

qc0A_C0= Qv0A_C0.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;
qc0N_C0= Qv0N_C0.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;
qc0N_C0_57= Qv0N_C0_57.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;
qc0N_C1= Qv0N_C1.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;

format short

 NIST = readmatrix('NISTSaturationProperties_Water.csv');
  c_l = interp1(NIST(:,1),NIST(:,10),T)
  c_v = interp1(NIST(:,1),NIST(:,22),T)

  T
  grooveDepth = Rg_star*1000
Re_v_CaseA  = max(2.*qc_CaseA./(pi*mu_v*RL_star*h_lv),[],2)'
Re_v_CaseB  = max(2.*qc_CaseB./(pi*mu_v*RL_star*h_lv),[],2)'
Re_v_CaseC0 = max(2.*qc0A_C0./(pi*mu_v*RL_star*h_lv),[],2)'
Re_v_CaseC1 = max(2.*qc0N_C1./(pi*mu_v*RL_star*h_lv),[],2)'
Ma_v_CaseA  = mu_v * Re_v_CaseA / (rho_v * RL_star * c_v)
Ma_v_CaseB  = mu_v * Re_v_CaseB / (rho_v * RL_star * c_v)
Ma_v_CaseC0 = mu_v * Re_v_CaseC0 / (rho_v * RL_star * c_v)
Ma_v_CaseC1 = mu_v * Re_v_CaseC1 / (rho_v * RL_star * c_v)

Re_l_CaseA  = 2.*max(qc_CaseA ,[],2)'.*tp./(pi*mu_l*h_lv .* (2*(Rg_star-RL_star) + tg*(Rg_star+RL_star)) )
Re_l_CaseB  = 2.*max(qc_CaseB ,[],2)'.*tp./(pi*mu_l*h_lv .* (2*(Rg_star-RL_star) + tg*(Rg_star+RL_star)) )
Re_l_CaseC0 = 2.*max(qc0A_C0  ,[],2)'.*tp./(pi*mu_l*h_lv .* (2*(Rg_star-RL_star) + tg*(Rg_star+RL_star)) )
Re_l_CaseC1 = 2.*max(qc0N_C1  ,[],2)'.*tp./(pi*mu_l*h_lv .* (2*(Rg_star-RL_star) + tg*(Rg_star+RL_star)) )
Ma_l_CaseA  = mu_l * Re_l_CaseA ./ (rho_l * Dh * c_l)
Ma_l_CaseB  = mu_l * Re_l_CaseB ./ (rho_l * Dh * c_l)
Ma_l_CaseC0 = mu_l * Re_l_CaseC0 ./ (rho_l * Dh * c_l)
Ma_l_CaseC1 = mu_l * Re_l_CaseC1 ./ (rho_l * Dh * c_l)


% 
plot(La_star,ans(3,:),'.');
title('200°C, Rg*=1.22mm')
xlabel('La* (m)')
ylabel('Re_v Case B');
box on; grid on;

% Save files
timestr = datetime('2025-07-17','Format','yyyy-MM-dd');
filename = sprintf('CapLimRgs%1.2i-%1.2i_T%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,timestr)
save(filename) % Saving Xi Numerical Results


ratioAC = (qc_CaseA(:,60) - qc0N_C1(:,60))./(qc0N_C1(:,60))

ratioC1A_range = (qc_CaseA-qc0N_C1)./(qc0N_C1);
ratioC1A_range_maxmin = [ min(ratioC1A_range,[],2) max(ratioC1A_range,[],2)] * 100
ratioAC0_range = (qc0A_C0-qc_CaseA)./(qc_CaseA);
ratioAC0_range_maxmin = [ min(ratioAC0_range,[],2) max(ratioAC0_range,[],2)] * 100
ratioAC1_range = (qc0N_C1-qc_CaseA)./(qc_CaseA);
ratioAC1_range_maxmin = [ min(ratioAC1_range,[],2) max(ratioAC1_range,[],2)] * 100

qcLa015 = qc0N_C1(:,60)

% Plotting Section (Capillary Limit VERISION 1) -------------------------------------


CArr = ["#0072BD","#D95319","#77AC30"];

close all; 
hfig = figure; hold on; box on

% for ig = 1:N_geom
%     nameStr = 'Analytical $R_g^*$=' + compose("%.2f",Rg_star(ig)*1E3) + 'mm';
%     % nameStr = 'Analytical $R_g^*$=' + compose("%.2f",Rg_star(ig)*1E3) + 'mm';
%     plot(La_star,qc0A_C0(ig,:),'-',Color=CArr(ig),DisplayName=nameStr)
% end
% for ig = 1:N_geom
%     nameStr = 'Numerical $R_g^*$=' + compose("%.2f",Rg_star(ig)*1E3) + 'mm';
%     plot(La_star,qc0N_C0(ig,:),'.' ,Color=CArr(ig),DisplayName=nameStr)
% end
% for ig = 1:N_geom
%     nameStr = 'Numerical C1 $R_g^*$=' + compose("%.2f",Rg_star(ig)*1E3) + 'mm';
%     plot(La_star,qc0N_C1(ig,:),'o' ,Color=CArr(ig),DisplayName=nameStr)
% end

idx_la = 10:110;

lw = 1.2;

for ig = 1:1

    % nameStr = 'Numerical';
    % plot(La_star(idx_la),qc0N_C0(ig,idx_la),'.' ,Color=CArr(ig),DisplayName=nameStr)
    nameStr = 'A'; plot(nan,nan,'--' ,LineWidth=lw,Color='k',DisplayName=nameStr) % nan is dummy plot for making legend
    plot(La_star(idx_la),qc_CaseA(ig,idx_la),'--' ,LineWidth=lw,Color=CArr(ig),HandleVisibility='off')
    nameStr = 'B';    plot(nan,nan,'-.',LineWidth=lw ,Color='k',DisplayName=nameStr)
    plot(La_star(idx_la),qc_CaseB(ig,idx_la),'-.' ,LineWidth=lw,Color=CArr(ig),HandleVisibility='off')
    nameStr = 'C$_0$';    plot(nan,nan,'-' ,LineWidth=lw,Color='k',DisplayName=nameStr)
    plot(La_star(idx_la),qc0A_C0(ig,idx_la),'-',LineWidth=lw,Color=CArr(ig),HandleVisibility='off')
    nameStr = 'C';    plot(nan,nan,'.' ,LineWidth=lw,Color='k',DisplayName=nameStr)
    plot(La_star(idx_la),qc0N_C1(ig,idx_la),'.',LineWidth=lw,Color=CArr(ig),HandleVisibility='off')

end

for ig = 2:3
    plot(La_star(idx_la),qc_CaseA(ig,idx_la),'--' ,LineWidth=lw,Color=CArr(ig),HandleVisibility='off')
    plot(La_star(idx_la),qc_CaseB(ig,idx_la),'-.' ,LineWidth=lw,Color=CArr(ig),HandleVisibility='off')
    plot(La_star(idx_la),qc0A_C0(ig,idx_la),'-',LineWidth=lw,Color=CArr(ig),HandleVisibility='off')
    % plot(La_star(idx_la),qc0N_C0(ig,idx_la),'.' ,Color=CArr(ig),HandleVisibility='off')
    plot(La_star(idx_la),qc0N_C1(ig,idx_la),'.' ,LineWidth=lw,Color=CArr(ig),HandleVisibility='off')

end
% xlim([0.025 0.275])
% ylim([2 16])

switch T
    case  60; 
    case 200; ylim([2 22]); %0.005;%0.0009;%[0.01 1]; % 200C-0.1
    case 300; ylim([0 8])
end
xlabel('$L_{\mathrm{a}}^*$ (m)') % (Lt* = ',num2str(Lt_star),' m)'])
ylabel('$q^*_{\mathrm{cap}}$ (W)')
% title(['Capillary limit vs Adiabatic sec length at T = ',num2str(T),'°C'])
legend(Location="northeast",Box="off",Orientation="vertical",Interpreter='latex')

hfigname = ['Figures\CapLim' sprintf('_T%i',T)];

picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document
    % set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
print(hfig,hfigname,'-dpng','-painters')

% % Volumetric Flowrate
% CArr = ["#0072BD","#D95319","#EDB120"];
% hfig = figure(2); hold on; grid on; box on
% hfig.Position = [210 120 900 900*0.75];
% for ig = 1:N_geom
%     nameStr = 'Andrew C0 Rg = ' + compose("%1.2i",Rg_star(ig));
%     plot(La_star,Qv0A_C0(ig,:),'-',Color=CArr(ig),DisplayName=nameStr)
% end
% for ig = 1:N_geom
%     nameStr = 'Numerical C0 Rg = ' + compose("%1.2i",Rg_star(ig));
%     plot(La_star,Qv0N_C0(ig,:),'.' ,Color=CArr(ig),DisplayName=nameStr)
% end
% for ig = 1:N_geom
%     nameStr = 'Numerical C1 Rg = ' + compose("%1.2i",Rg_star(ig));
%     plot(La_star,Qv0N_C1(ig,:),'o' ,Color=CArr(ig),DisplayName=nameStr)
% end
% xlabel(['La* [m] (Lt* = ',num2str(Lt_star),' m)'])
% ylabel('Qv0 [W]')
% title(['Volumetric flowrate vs Adiabatic sec length at T = ',num2str(T),'°C'])
% legend(Location="northwest")


%% p* vs z* Deep Dive Study 
clear all; close all; clc; format  compact; format long
load('CapLimRgs1.18e+00-1.22e+00_T60_2025-07-17.mat')
% load('CapLimRgs1.18e+00-1.22e+00_T200_2025-07-17.mat')
% load('CapLimRgs1.18e+00-1.22e+00_T300_2025-07-17.mat')

N_pt = 1000;

ig = 1; % Grooved depth value
iz = 80; % La* location 40=0.2, 80=0.1, 100=0.05
La_star(iz)
i  = 15; % \theta_m^0 location
i_57 = 177 ; % \theta_m^0 = 57^o location

dxidz_v_N = Qv0N_C0(ig,iz)./tQv0N(ig,:);
dxidz_l_N = dxidz_v_N .* tXi0N(ig,:);
dxidz_v_A = Qv0A_C0(ig,iz)./tQv0A(ig,1);
dxidz_l_A = dxidz_v_A .* tXi0A(ig,1);

dxidz_v_N_57 = Qv0N_C0_57(ig,iz)./tQv0N(ig,:);
dxidz_l_N_57 = dxidz_v_N_57 .* tXi0N(ig,:);

ze_star = linspace(-Le_star(iz),                      0, 50);
zc_star = linspace( La_star(iz),La_star(iz)+Lc_star(iz), 50);
za      = linspace(           0,                      1, 50);
za_star = za.*La_star(iz);
z_star = [ze_star za_star zc_star];


for i=i % Numerical  Case C0
p_veN = 1 - 8*Ng*Qv0N_C0(ig,iz)*Le(iz)/(pi*RL^4) ...
         - 1*Ng*rho*Qv0N_C0(ig,iz)*Le(iz)/(mu*Aw(ig))+d_Bo(iz)*Le(iz);
d_pcN = 8*Ng*Qv0N_C0(ig,iz)*Lc(iz)/(pi*RL^4) ...
         + 1*Ng*rho*Qv0N_C0(ig,iz)*Lc(iz)/(mu*Aw(ig))-d_Bo(iz)*Lc(iz);

d_pe_starN =  p_veN .* Dp_cmax_star;
d_pc_starN =  d_pcN .* Dp_cmax_star;

% Evaporator-----------
p_ve_starN = p_sat - 8*mu_v*qc0N_C0(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) ...
    - rho_v*g*sin(phi)*(Le_star(iz)+ze_star);
p_le_starN = p_sat - Dp_cmax_star + mu_l*qc0N_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2)...
    - rho_l*g*sin(phi)*(Le_star(iz)+ze_star);
diffcheck_evap = p_ve_starN(end) - p_le_starN(end) - d_pe_starN; 

% Adiabatic  -----------
p_vaN = (-Bo_v(iz) - dxidz_v_N(i)) .* za + d_pe_starN/Dp_cmax_star; % pva
p_laN = (-Bo_l(iz) - dxidz_l_N(i)) .* za;
p_va_starN = p_vaN.*Dp_cmax_star + p_le_starN(end);
p_la_starN = p_laN.*Dp_cmax_star + p_le_starN(end);
diffcheck_adia = p_va_starN(end) - p_la_starN(end) - d_pc_starN ;

% Condnesor  -----------
p_vc_starN = p_va_starN(end) - 8*mu_v*qc0N_C0(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_v*g*sin(phi)*(zc_star-La_star(iz));
p_lc_starN = p_la_starN(end) + mu_l*qc0N_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_l*g*sin(phi)*(zc_star-La_star(iz));

p_v_starN = [p_ve_starN p_va_starN p_vc_starN];
p_l_starN = [p_le_starN p_la_starN p_lc_starN];
end 

for i_57=i_57 % Numerical  Case C0_57
p_veN_57 = 1 - 8*Ng*Qv0N_C0_57(ig,iz)*Le(iz)/(pi*RL^4) ...
         - 1*Ng*rho*Qv0N_C0_57(ig,iz)*Le(iz)/(mu*Aw(ig))+d_Bo(iz)*Le(iz);
d_pcN_57 = 8*Ng*Qv0N_C0_57(ig,iz)*Lc(iz)/(pi*RL^4) ...
         + 1*Ng*rho*Qv0N_C0_57(ig,iz)*Lc(iz)/(mu*Aw(ig))-d_Bo(iz)*Lc(iz);

d_pe_starN_57 =  p_veN_57 .* Dp_cmax_star;
d_pc_starN_57 =  d_pcN_57 .* Dp_cmax_star;

% Evaporator-----------
p_ve_starN_57 = p_sat - 8*mu_v*qc0N_C0_57(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) ...
    - rho_v*g*sin(phi)*(Le_star(iz)+ze_star);
p_le_starN_57 = p_sat - Dp_cmax_star + mu_l*qc0N_C0_57(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2)...
    - rho_l*g*sin(phi)*(Le_star(iz)+ze_star);
diffcheck_evap = p_ve_starN_57(end) - p_le_starN_57(end) - d_pe_starN_57; 

% Adiabatic  -----------
p_vaN_57 = (-Bo_v(iz) - dxidz_v_N_57(i_57)) .* za + d_pe_starN_57/Dp_cmax_star; % pva
p_laN_57 = (-Bo_l(iz) - dxidz_l_N_57(i_57)) .* za;
p_va_starN_57 = p_vaN_57.*Dp_cmax_star + p_le_starN_57(end);
p_la_starN_57 = p_laN_57.*Dp_cmax_star + p_le_starN_57(end);
diffcheck_adia_57 = p_va_starN_57(end) - p_la_starN_57(end) - d_pc_starN_57;

% Condnesor  -----------
p_vc_starN_57  = p_va_starN_57 (end) - 8*mu_v*qc0N_C0_57(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_v*g*sin(phi)*(zc_star-La_star(iz));
p_lc_starN_57  = p_la_starN_57 (end) + mu_l*qc0N_C0_57 (ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_l*g*sin(phi)*(zc_star-La_star(iz));

p_v_starN_57 = [p_ve_starN_57 p_va_starN_57 p_vc_starN_57];
p_l_starN_57 = [p_le_starN_57 p_la_starN_57 p_lc_starN_57];

dpvdz_starN_Adia_57 = (p_va_starN_57(1) - p_va_starN_57(end))/La_star(iz);
dpldz_starN_Adia_57 = (p_la_starN_57(1) - p_la_starN_57(end))/La_star(iz);

tot_p_grad_CaseC0_57 = abs(dpvdz_starN_Adia_57) + abs(dpldz_starN_Adia_57)
end 

for i=i % Analytical Case C0
p_veA = 1 - 8*Ng*Qv0A_C0(ig,iz)*Le(iz)/(pi*RL^4) ...
         - 1*Ng*rho*Qv0A_C0(ig,iz)*Le(iz)/(mu*Aw(ig))+d_Bo(iz)*Le(iz);
d_pcA = 8*Ng*Qv0A_C0(ig,iz)*Lc(iz)/(pi*RL^4) ...
         + 1*Ng*rho*Qv0A_C0(ig,iz)*Lc(iz)/(mu*Aw(ig))-d_Bo(iz)*Lc(iz);

d_pe_starA =  p_veA .* Dp_cmax_star;
d_pc_starA =  d_pcA .* Dp_cmax_star;


% Evaporator-----------
p_ve_starA = p_sat - 8*mu_v*qc0A_C0(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) ...
    - rho_v*g*sin(phi)*(Le_star(iz)+ze_star);
p_le_starA = p_sat - Dp_cmax_star + mu_l*qc0A_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2)...
    - rho_l*g*sin(phi)*(Le_star(iz)+ze_star);
diffcheck_evap = p_ve_starA(end) - p_le_starA(end) - d_pe_starA ;

% Adiabatic  -----------
p_vaA = (-Bo_v(iz) + dxidz_v_A) .* za + d_pe_starA/Dp_cmax_star; % pva
p_laA = (-Bo_l(iz) + dxidz_l_A) .* za;
p_va_starA = p_vaA.*Dp_cmax_star + p_le_starA(end);
p_la_starA = p_laA.*Dp_cmax_star + p_le_starA(end);
diffcheck_adia = p_va_starA(end) - p_la_starA(end) - d_pc_starA; 

% Condnesor  -----------
p_vc_starA = p_va_starA(end) - 8*mu_v*qc0A_C0(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_v*g*sin(phi)*(zc_star-La_star(iz));
p_lc_starA = p_la_starA(end) + mu_l*qc0A_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_l*g*sin(phi)*(zc_star-La_star(iz));

p_v_starA = [p_ve_starA p_va_starA p_vc_starA];
p_l_starA = [p_le_starA p_la_starA p_lc_starA];

dpvdz_starA_Adia = (p_va_starA(1) - p_va_starA(end))/La_star(iz);
dpldz_starA_Adia = (p_la_starA(1) - p_la_starA(end))/La_star(iz);

tot_p_grad_CaseC0 = abs(dpvdz_starA_Adia) + abs(dpldz_starA_Adia)

end 

for i=i  % Case A & C0 comparision

delta_p_vt_AC = 8*mu_v*qc0A_C0(ig,iz)*Leff_star(iz)/(pi*rho_v*RL_star^4*h_lv)+...
    rho_v*g*sin(phi)*Lt_star;

p_ve_star_CaseAC = p_sat - 8*mu_v*qc0A_C0(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) - rho_v*g*sin(phi)*...
    (Le_star(iz)+ze_star);
p_va_star_CaseAC = p_sat - 8*mu_v*qc0A_C0(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (za_star+Le_star(iz)/2) - rho_v*g*sin(phi)*(Le_star(iz)+za_star);
p_vc_star_CaseAC = p_sat - 8*mu_v*qc0A_C0(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (Le_star(iz)/2+La_star(iz)+(zc_star-La_star(iz)).*(La_star(iz)+2*Lc_star(iz)-zc_star)/...
    (2*Lc_star(iz))) - rho_v*g*sin(phi)*(Le_star(iz)+zc_star);

p_le_star_CaseAC = p_sat - delta_p_vt_AC + ...
    mu_l*qc0A_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz)) + ze_star - Lc_star(iz)/2 - La_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - ze_star);
p_la_star_CaseAC = p_sat - delta_p_vt_AC + ...
    mu_l*qc0A_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (za_star - Lc_star(iz)/2 - La_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - za_star);
p_lc_star_CaseAC = p_sat - delta_p_vt_AC - ...
    mu_l*qc0A_C0(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (La_star(iz) + Lc_star(iz) - zc_star).^2/(2*Lc_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - zc_star);

z_star = [ze_star za_star zc_star];
p_v_star_CaseAC = [p_ve_star_CaseAC p_va_star_CaseAC p_vc_star_CaseAC];
p_l_star_CaseAC = [p_le_star_CaseAC p_la_star_CaseAC p_lc_star_CaseAC];

dpvdz_star_CaseAC_Adia = (p_va_star_CaseAC(1) - p_va_star_CaseAC(end))/La_star(iz);
dpldz_star_CaseAC_Adia = (p_la_star_CaseAC(1) - p_la_star_CaseAC(end))/La_star(iz);

tot_p_grad_CaseAC = abs(dpvdz_star_CaseAC_Adia) + abs(dpldz_star_CaseAC_Adia)

ratio_caseAC = tot_p_grad_CaseAC/tot_p_grad_CaseC0

end 

for i=i  % Case A & C0_57 comparision

delta_p_vt_AC_57 = 8*mu_v*qc0N_C0_57(ig,iz)*Leff_star(iz)/(pi*rho_v*RL_star^4*h_lv)+...
    rho_v*g*sin(phi)*Lt_star;

p_ve_star_CaseAC_57 = p_sat - 8*mu_v*qc0N_C0_57(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) - rho_v*g*sin(phi)*...
    (Le_star(iz)+ze_star);
p_va_star_CaseAC_57 = p_sat - 8*mu_v*qc0N_C0_57(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (za_star+Le_star(iz)/2) - rho_v*g*sin(phi)*(Le_star(iz)+za_star);
p_vc_star_CaseAC_57 = p_sat - 8*mu_v*qc0N_C0_57(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (Le_star(iz)/2+La_star(iz)+(zc_star-La_star(iz)).*(La_star(iz)+2*Lc_star(iz)-zc_star)/...
    (2*Lc_star(iz))) - rho_v*g*sin(phi)*(Le_star(iz)+zc_star);

p_le_star_CaseAC_57 = p_sat - delta_p_vt_AC + ...
    mu_l*qc0N_C0_57(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz)) + ze_star - Lc_star(iz)/2 - La_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - ze_star);
p_la_star_CaseAC_57 = p_sat - delta_p_vt_AC + ...
    mu_l*qc0N_C0_57(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (za_star - Lc_star(iz)/2 - La_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - za_star);
p_lc_star_CaseAC_57 = p_sat - delta_p_vt_AC - ...
    mu_l*qc0N_C0_57(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (La_star(iz) + Lc_star(iz) - zc_star).^2/(2*Lc_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - zc_star);

z_star = [ze_star za_star zc_star];
p_v_star_CaseAC_57 = [p_ve_star_CaseAC_57 p_va_star_CaseAC_57 p_vc_star_CaseAC_57];
p_l_star_CaseAC_57 = [p_le_star_CaseAC_57 p_la_star_CaseAC_57 p_lc_star_CaseAC_57];

dpvdz_star_CaseAC_Adia_57 = (p_va_star_CaseAC_57(1) - p_va_star_CaseAC_57(end))/La_star(iz);
dpldz_star_CaseAC_Adia_57 = (p_la_star_CaseAC_57(1) - p_la_star_CaseAC_57(end))/La_star(iz);

tot_p_grad_CaseAC_57 = abs(dpvdz_star_CaseAC_Adia_57) + abs(dpldz_star_CaseAC_Adia_57)

ratio_caseAC_57 = tot_p_grad_CaseAC_57/tot_p_grad_CaseC0_57

end 

for i=i % Numerical Case C1
    
p_veN1 = 1 - 8*Ng*Qv0N_C1(ig,iz)*Le(iz)/(pi*RL^4) ...
         - 1*Ng*rho*Qv0N_C1(ig,iz)*Le(iz)/(mu*Aw(ig))+d_Bo(iz)*Le(iz);
d_pcN1 = 8*Ng*Qv0N_C1(ig,iz)*Lc(iz)/(pi*RL^4) ...
         + 1*Ng*rho*Qv0N_C1(ig,iz)*Lc(iz)/(mu*Aw(ig))-d_Bo(iz)*Lc(iz);

d_pe_starN1 =  p_veN1 .* Dp_cmax_star;
d_pc_starN1 =  d_pcN1 .* Dp_cmax_star;

% Evaporator-----------
p_ve_starN1 = p_sat - 8*mu_v*qc0N_C1(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) ...
    - rho_v*g*sin(phi)*(Le_star(iz)+ze_star);
p_le_starN1 = p_sat - Dp_cmax_star + mu_l*qc0N_C1(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2)...
    - rho_l*g*sin(phi)*(Le_star(iz)+ze_star);

diffcheck_evap_N1 = p_ve_starN1(end) - p_le_starN1(end) - d_pe_starN1 ;

% Adiabatic  -----------
% p_vaN1 = (-Bo_v(iz) - dxidz_v_N(i)) .* za + d_pe_starN1/Dp_cmax_star; % pva
% p_laN1 = (-Bo_l(iz) - dxidz_l_N(i)) .* za;

% Case C1 

Rm0_e = 1/(Ca_hat*p_veN1);
Rm0_c = 1/(Ca_hat*d_pcN1);
tm0_e = interp1(Rm0,tm0,Rm0_e)
tm0_c = interp1(Rm0,tm0,Rm0_c) 

tm0q = linspace(tm0_e,tm0_c,N_pt);
Rm0q = interp1(tm0,Rm0,tm0q);
delta_p0 = 1./(Rm0q*Ca_hat);

tXi0Nq = interp1(tm0,tXi0N(ig,:),tm0q);
tQv0Nq = interp1(tm0,tQv0N(ig,:),tm0q);

dxidz_v_Nq = interp1(tm0,dxidz_v_N,tm0q);
dxidz_l_Nq = dxidz_v_Nq .* tXi0Nq;
zaq = zeros(size(delta_p0));
dp_starN1 = delta_p0.*Dp_cmax_star;

zaq = cumtrapz(delta_p0,...
                1./(d_Bo(iz)+Qv0N_C1(ig,iz).*(1-tXi0Nq)./(-tQv0Nq)));
zaq_star = zaq.*La_star(iz);

p_vaNq = cumtrapz(zaq, -Bo_v(iz)- dxidz_v_Nq) + d_pe_starN1/Dp_cmax_star;
% p_laNq = cumtrapz(zaq, -Bo_l(iz)- dxidz_l_Nq);
p_va_starNq = p_vaNq.*Dp_cmax_star + p_le_starN1(end);
% p_la_starNq = p_laNq.*Dp_cmax_star + p_le_starN1(end);

% p_va_starNq = p_la_starNq + dp_starN1;
p_la_starNq = p_va_starNq - dp_starN1;

diffcheck_adiaN1 = p_va_starNq(end) - p_la_starNq(end) - d_pc_starN1;

% Condnesor  -----------
p_vc_starN1 = p_va_starNq(end) - 8*mu_v*qc0N_C1(ig,iz)./(pi*rho_v*RL_star^4*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_v*g*sin(phi)*(zc_star-La_star(iz));
p_lc_starN1 = p_la_starNq(end) + mu_l*qc0N_C1(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
      ( (La_star(iz)+Lc_star(iz)).*zc_star - zc_star.^2/2 - La_star(iz).^2/2 - La_star(iz)*Lc_star(iz))/Lc_star(iz)... 
      - rho_l*g*sin(phi)*(zc_star-La_star(iz));


p_va_starN = interp1(zaq_star,p_va_starNq,za_star);
p_vl_starN = interp1(zaq_star,p_la_starNq,za_star);

z_starN1   = [ze_star za_star zc_star];
p_v_starN1 = [p_ve_starN1 p_va_starN p_vc_starN1];
p_l_starN1 = [p_le_starN1 p_vl_starN p_lc_starN1];

end 	

for i=i  % Analytical  Case A
delta_p_vt_A = 8*mu_v*qc_CaseA(ig,iz)*Leff_star(iz)/(pi*rho_v*RL_star^4*h_lv)+...
    rho_v*g*sin(phi)*Lt_star;

p_ve_star_A = p_sat - 8*mu_v*qc_CaseA(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) - rho_v*g*sin(phi)*...
    (Le_star(iz)+ze_star);
p_va_star_A = p_sat - 8*mu_v*qc_CaseA(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (za_star+Le_star(iz)/2) - rho_v*g*sin(phi)*(Le_star(iz)+za_star);
p_vc_star_A = p_sat - 8*mu_v*qc_CaseA(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (Le_star(iz)/2+La_star(iz)+(zc_star-La_star(iz)).*(La_star(iz)+2*Lc_star(iz)-zc_star)/...
    (2*Lc_star(iz))) - rho_v*g*sin(phi)*(Le_star(iz)+zc_star);

p_le_star_A = p_sat - delta_p_vt_A + ...
    mu_l*qc_CaseA(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star(iz)) + ze_star - Lc_star(iz)/2 - La_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - ze_star);
p_la_star_A = p_sat - delta_p_vt_A + ...
    mu_l*qc_CaseA(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (za_star - Lc_star(iz)/2 - La_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - za_star);
p_lc_star_A = p_sat - delta_p_vt_A - ...
    mu_l*qc_CaseA(ig,iz)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (La_star(iz) + Lc_star(iz) - zc_star).^2/(2*Lc_star(iz)) + ...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - zc_star);

z_star = [ze_star za_star zc_star];
p_v_star_CaseA = [p_ve_star_A p_va_star_A p_vc_star_A];
p_l_star_CaseA = [p_le_star_A p_la_star_A p_lc_star_A];
end 

for i=i  % Analytical  Case B
[K1_star,K2_star,K4_star,Ca_star] = A_FullMassFlowRate_CaseAB(mu,tp,tg,RL_star,RL,Rg_star(ig),RgArr(ig),P_star,tXi0A(ig,1),tXi0A(ig,2),N_An);

delta_p_vt_B = 8*mu_v*qc_CaseB(ig,iz)*Leff_star(iz)/(pi*rho_v*RL_star^4*h_lv)+...
    rho_v*g*sin(phi)*Lt_star; %is this right???

p_ve_star_B = p_sat - 8*mu_v*qc_CaseB(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star(iz))+ze_star+Le_star(iz)/2) - rho_v*g*sin(phi)*...
    (Le_star(iz)+ze_star);
p_va_star_B = p_sat - 8*mu_v*qc_CaseB(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (za_star+Le_star(iz)/2) - rho_v*g*sin(phi)*(Le_star(iz)+za_star);
p_vc_star_B = p_sat - 8*mu_v*qc_CaseB(ig,iz)/(pi*rho_v*RL_star^4*h_lv)*...
    (Le_star(iz)/2+La_star(iz)+(zc_star-La_star(iz)).*(La_star(iz)+2*Lc_star(iz)-zc_star)/...
    (2*Lc_star(iz))) - rho_v*g*sin(phi)*(Le_star(iz)+zc_star);

p_le_star_B = p_sat - delta_p_vt_B + ...
    (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv))...
    *mu_l*qc_CaseB(ig,iz)/(K1_star + K2_star - Ca_star*K4_star).*...
    (ze_star.^2./(2*Le_star(iz)) + ze_star - Lc_star(iz)/2 - La_star(iz))+...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - ze_star);
p_la_star_B = p_sat - delta_p_vt_B + ...
    (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv))...
    *mu_l*qc_CaseB(ig,iz)/(K1_star + K2_star - Ca_star*K4_star)*...
    (za_star - Lc_star(iz)/2 - La_star(iz))+...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - za_star);
p_lc_star_B = p_sat - delta_p_vt_B - ...
    (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv))...
    *mu_l*qc_CaseB(ig,iz)/(K1_star + K2_star - Ca_star*K4_star)*...
    (La_star(iz) + Lc_star(iz) - zc_star).^2/(2*Lc_star(iz))+...
    rho_l*g*sin(phi)*(La_star(iz) + Lc_star(iz) - zc_star);

p_v_star_CaseB = [p_ve_star_B p_va_star_B p_vc_star_B];
p_l_star_CaseB = [p_le_star_B p_la_star_B p_lc_star_B];

end

timestr = datetime('2025-11-27','Format','yyyy-MM-dd');
filename = sprintf('pzplot_T%i_Rgs%1.2i_La_star%1.2i_%s.mat',T,Rg_star(ig)*1E3,La_star(iz),timestr)
save(['Pressure_Plots\' filename]) % Saving Xi Numerical Results

hfig = figure; % save the figure handle in a variable
hold on; grid on; box on
plot(z_starN1,p_v_starN1*0.001,'.',Color='#0072BD',MarkerSize=3,DisplayName='$p^{*}$');
plot(z_starN1,p_l_starN1*0.001,'.',Color='#0072BD',MarkerSize=3,HandleVisibility='off');
plot(z_star  ,p_v_starA*0.001,'-',Color='#D95319'    ,LineWidth=1,DisplayName='$p^{*^{0}}$');
plot(z_star  ,p_l_starA*0.001,'-',Color='#D95319'    ,LineWidth=1,HandleVisibility='off');
% plot(z_star  ,p_v_starN*0.001,'.',Color='#D95319',MarkerSize=6,DisplayName='Numerical');
% plot(z_star  ,p_l_starN*0.001,'.',Color='#D95319',MarkerSize=6,HandleVisibility='off');
plot(z_starN1,p_v_star_CaseA*0.001,'--',Color='#EDB120',MarkerSize=3,DisplayName='$p^{*^{\rm{A}}}$');
plot(z_starN1,p_l_star_CaseA*0.001,'--',Color='#EDB120',MarkerSize=3,HandleVisibility='off');
plot(z_starN1,p_v_star_CaseB*0.001,'-.',Color='#7E2F8E',MarkerSize=3,DisplayName='$p^{*^{\rm{B}}}$');
plot(z_starN1,p_l_star_CaseB*0.001,'-.',Color='#7E2F8E',MarkerSize=3,HandleVisibility='off');
% xline(0,'-',{'Adiabatic','begin'},LabelVerticalAlignment='bottom');
% xline(La_star(iz),HandleVisibility='off')
xlabel('$z^*$ (m)'); ylabel('$p^*$ (kPa)'); 
% xlim([-0.1 -0.1+Lt_star])
legend(Location='southeast')
title([ 'T = ',num2str(T),'$^\circ$C;    ',...
        'Rg* = ',num2str(Rg_star(ig)*1000),'mm;    ',...
        'La* = ',num2str(La_star(iz)),'m;    ',...
        'ACratio* = ',num2str(ratio_caseAC),';    '])

    % txt = ['$toto presstio gradient ratio~=~$' sprintf('%1.4f',ratio_caseAC)];
    % text(min(z_star),min(p_l_starA),txt,'FontSize',12,Interpreter='latex')

hfigname = ['Figures\p_vs_z' sprintf('_T%i',T) '_Rgstr' num2str(Rg_star(ig)*100000) '_Lastr' num2str(La_star(iz)*100)];

picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 1; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
% print(hfig,hfigname,'-dpng','-painters')


%% Plotting Section (p* vs z* VERISION 2) -------------------------------------
clear all; close all; clc; format  compact; format long

filenameCell{1,1} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
filenameCell{1,2} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 
filenameCell{2,1} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
filenameCell{2,2} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 
filenameCell{3,1} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
filenameCell{3,2} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 
titleCell{1,1} ='$T_{\rm{o}}^*=60^{\circ}$C, $L_{\rm{a}}^*=0.1$~m';
titleCell{1,2} ='$T_{\rm{o}}^*=60^{\circ}$C, $L_{\rm{a}}^*=0.2$~m';
titleCell{2,1} ='$T_{\rm{o}}^*=200^{\circ}$C, $L_{\rm{a}}^*=0.1$~m';
titleCell{2,2} ='$T_{\rm{o}}^*=200^{\circ}$C, $L_{\rm{a}}^*=0.2$~m';
titleCell{3,1} ='$T_{\rm{o}}^*=300^{\circ}$C, $L_{\rm{a}}^*=0.1$~m';
titleCell{3,2} ='$T_{\rm{o}}^*=300^{\circ}$C, $L_{\rm{a}}^*=0.2$~m';

hfig = figure('DefaultAxesFontSize',12); % save the figure handle in a variable
fsLgd = 12*1.5
fsAxi = 14*1.5

ls = 1.2
i_sub = 0

tcl = tiledlayout(3,2,"Padding","tight")


for i_colum = 1:3
for i_row = 1:2

    i_sub = i_sub + 1;
    load(filenameCell{i_colum,i_row})
    nexttile
    hold on; grid on; box on
    plot(z_starN1,p_v_starN1*1E-6,'.',Color='#0072BD',DisplayName='$p^{*}$');
    plot(z_starN1,p_l_starN1*1E-6,'.',Color='#0072BD',HandleVisibility='off');
    plot(z_star  ,p_v_starA*1E-6,'-',Color='#D95319',LineWidth=ls,DisplayName='$p^{*^{0}}$');
    plot(z_star  ,p_l_starA*1E-6,'-',Color='#D95319',LineWidth=ls,HandleVisibility='off');
    % plot(z_star  ,p_v_starN*1E-6,'.',Color='#D95319',MarkerSize=6,DisplayName='Numerical');
    % plot(z_star  ,p_l_starN*1E-6,'.',Color='#D95319',MarkerSize=6,HandleVisibility='off');
    plot(z_starN1,p_v_star_CaseA*1E-6,'--',Color='#77AC30',LineWidth=ls,DisplayName='$p^{*^{\rm{A}}}$');
    plot(z_starN1,p_l_star_CaseA*1E-6,'--',Color='#77AC30',LineWidth=ls,HandleVisibility='off');
    plot(z_starN1,p_v_star_CaseB*1E-6,'-.',Color='#EDB120',LineWidth=ls,DisplayName='$p^{*^{\rm{B}}}$');
    plot(z_starN1,p_l_star_CaseB*1E-6,'-.',Color='#EDB120',LineWidth=ls,HandleVisibility='off')
    lgd = legend(Location='best',Orientation='horizontal',Box='off',FontSize=fsLgd)    
    title(lgd,titleCell{i_colum,i_row})
    titleCell{i_colum,i_row}
    capPressure = p_v_starN1(50) - p_l_starN1(50)

    if i_colum == 3; xlabel('$z^*$ (m)',FontSize=fsAxi);
    end 
    if i_row == 1; xlim([-0.1 0.2]); ylabel('$p^*$ (MPa)',FontSize=fsAxi); 
    else xlim([-0.05 0.25]);
    end
    
% xlim([-0.1 -0.1+Lt_star])

end 
end

hfigname = ['Figures\p_star_vs_z_star_v2'];
picturewidth = 35; % set this parameter and keep it forever
hw_ratio =0.6; % feel free to play with this ratio
set(findall(hfig,'-property','MarkerSize'),'MarkerSize',8)
% set(findall(hfig,'-property','LineWidth'),'LineWidth',1.2)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
movegui('south')
print(hfig,hfigname,'-dpdf','-painters')
print(hfig,hfigname,'-dpng','-painters')
%% Plotting Section (p* vs z* VERISION 3, Paper Plot) -------------------------------------
clear all; close all; clc; format  compact; format long

titleCell{1,1} ='$T_{\rm{o}}^*=60^{\circ}$C, $L_{\rm{a}}^*=0.1$~m';
titleCell{1,2} ='$T_{\rm{o}}^*=60^{\circ}$C, $L_{\rm{a}}^*=0.2$~m';
titleCell{2,1} ='$T_{\rm{o}}^*=200^{\circ}$C, $L_{\rm{a}}^*=0.1$~m';
titleCell{2,2} ='$T_{\rm{o}}^*=200^{\circ}$C, $L_{\rm{a}}^*=0.2$~m';
titleCell{3,1} ='$T_{\rm{o}}^*=300^{\circ}$C, $L_{\rm{a}}^*=0.1$~m';
titleCell{3,2} ='$T_{\rm{o}}^*=300^{\circ}$C, $L_{\rm{a}}^*=0.2$~m';

hfigname = ['Figures\p_star_vs_z_star_v3_Rg122'];
filenameCell{1,1} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
filenameCell{1,2} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 
filenameCell{2,1} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
filenameCell{2,2} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 
filenameCell{3,1} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
filenameCell{3,2} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 

% hfigname = ['Figures\p_star_vs_z_star_v3_Rg120'];
% filenameCell{1,1} = 'Pressure_Plots\pzplot_T60_Rgs1.20e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{1,2} = 'Pressure_Plots\pzplot_T60_Rgs1.20e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{2,1} = 'Pressure_Plots\pzplot_T200_Rgs1.20e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{2,2} = 'Pressure_Plots\pzplot_T200_Rgs1.20e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{3,1} = 'Pressure_Plots\pzplot_T300_Rgs1.20e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{3,2} = 'Pressure_Plots\pzplot_T300_Rgs1.20e+00_La_star2.00e-01_2025-07-17'; 

% hfigname = ['Figures\p_star_vs_z_star_v3_Rg118'];
% filenameCell{1,1} = 'Pressure_Plots\pzplot_T60_Rgs1.18e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{1,2} = 'Pressure_Plots\pzplot_T60_Rgs1.18e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{2,1} = 'Pressure_Plots\pzplot_T200_Rgs1.18e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{2,2} = 'Pressure_Plots\pzplot_T200_Rgs1.18e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{3,1} = 'Pressure_Plots\pzplot_T300_Rgs1.18e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{3,2} = 'Pressure_Plots\pzplot_T300_Rgs1.18e+00_La_star2.00e-01_2025-07-17'; 



% titleCell{1,1} ='$R_{\rm{g}}^*=1.18$~mm, $L_{\rm{a}}^*=0.1$~m';
% titleCell{1,2} ='$R_{\rm{g}}^*=1.18$~mm, $L_{\rm{a}}^*=0.2$~m';
% titleCell{2,1} ='$R_{\rm{g}}^*=1.20$~mm, $L_{\rm{a}}^*=0.1$~m';
% titleCell{2,2} ='$R_{\rm{g}}^*=1.20$~mm, $L_{\rm{a}}^*=0.2$~m';
% titleCell{3,1} ='$R_{\rm{g}}^*=1.22$~mm, $L_{\rm{a}}^*=0.1$~m';
% titleCell{3,2} ='$R_{\rm{g}}^*=1.22$~mm, $L_{\rm{a}}^*=0.2$~m';
% 
% hfigname = ['Figures\p_star_vs_z_star_v3_T60'];
% filenameCell{1,1} = 'Pressure_Plots\pzplot_T60_Rgs1.18e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{1,2} = 'Pressure_Plots\pzplot_T60_Rgs1.18e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{2,1} = 'Pressure_Plots\pzplot_T60_Rgs1.20e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{2,2} = 'Pressure_Plots\pzplot_T60_Rgs1.20e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{3,1} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{3,2} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 

% hfigname = ['Figures\p_star_vs_z_star_v3_T200'];
% filenameCell{1,1} = 'Pressure_Plots\pzplot_T200_Rgs1.18e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{1,2} = 'Pressure_Plots\pzplot_T200_Rgs1.18e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{2,1} = 'Pressure_Plots\pzplot_T200_Rgs1.20e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{2,2} = 'Pressure_Plots\pzplot_T200_Rgs1.20e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{3,1} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{3,2} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 

% hfigname = ['Figures\p_star_vs_z_star_v3_T300'];
% filenameCell{1,1} = 'Pressure_Plots\pzplot_T300_Rgs1.18e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{1,2} = 'Pressure_Plots\pzplot_T300_Rgs1.18e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{2,1} = 'Pressure_Plots\pzplot_T300_Rgs1.20e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{2,2} = 'Pressure_Plots\pzplot_T300_Rgs1.20e+00_La_star2.00e-01_2025-07-17'; 
% filenameCell{3,1} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star1.00e-01_2025-07-17'; 
% filenameCell{3,2} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star2.00e-01_2025-07-17'; 

hfig = figure('DefaultAxesFontSize',12); % save the figure handle in a variable
fsLgd = 12*1.5
fsAxi = 14*1.5

ls = 1.2
i_sub = 0

tcl = tiledlayout(4,2,"Padding","tight")


for i_colum = 1:3
for i_row = 1:2

    i_sub = i_sub + 1;
    load(filenameCell{i_colum,i_row})

    tm0_Cell{i_colum,i_row} = interp1(zaq_star,tm0q,za_star);
    z_star_Cell{i_colum,i_row} = za_star;

    nexttile
    hold on; grid on; box on
    plot(z_starN1,p_v_starN1*1E-6,'.',Color='#0072BD',DisplayName='C');
    plot(z_starN1,p_l_starN1*1E-6,'.',Color='#0072BD',HandleVisibility='off');
    plot(z_star  ,p_v_starA*1E-6,'-',Color='#D95319',LineWidth=ls,DisplayName='C$_0$');
    plot(z_star  ,p_l_starA*1E-6,'-',Color='#D95319',LineWidth=ls,HandleVisibility='off');
    % plot(z_star  ,p_v_starN*1E-6,'.',Color='#D95319',MarkerSize=6,DisplayName='Numerical');
    % plot(z_star  ,p_l_starN*1E-6,'.',Color='#D95319',MarkerSize=6,HandleVisibility='off');
    plot(z_starN1,p_v_star_CaseA*1E-6,'--',Color='#77AC30',LineWidth=ls,DisplayName='A');
    plot(z_starN1,p_l_star_CaseA*1E-6,'--',Color='#77AC30',LineWidth=ls,HandleVisibility='off');
    plot(z_starN1,p_v_star_CaseB*1E-6,'-.',Color='black',LineWidth=ls,DisplayName='B');
    plot(z_starN1,p_l_star_CaseB*1E-6,'-.',Color='black',LineWidth=ls,HandleVisibility='off')
    lgd = legend(Location='best',Orientation='horizontal',Box='off',FontSize=fsLgd)    
    title(lgd,titleCell{i_colum,i_row})
    titleCell{i_colum,i_row}
    capPressure = p_v_starN1(50) - p_l_starN1(50)

    % if i_colum == 3; xlabel('$z^*$ (m)',FontSize=fsAxi);
    % end 
    if i_row == 1; xlim([-0.1 0.2]); ylabel('$p^*$ (MPa)',FontSize=fsAxi); 
    else xlim([-0.05 0.25]);
    end
    


% xlim([-0.1 -0.1+Lt_star])

end 
end

    % tclText =  sprintf('$R_{\\rm{g}}^* = %1.2f\\rm{mm}$',Rg_star(ig)*1000)
    % tclText =  sprintf('$T_{\\rm{o}}^* = %i^\\circ$C',T)
    % title(tcl,tclText,FontSize = fsAxi,Interpreter='latex')
    nexttile
    hold on; grid on; box on
    plot(z_star_Cell{1,1},rad2deg(tm0_Cell{1,1}),'o',Color='black',DisplayName='$T_{\rm{o}}^*=60^{\circ}$C');
    plot(z_star_Cell{2,1},rad2deg(tm0_Cell{2,1}),'^',Color='black',DisplayName='$T_{\rm{o}}^*=200^{\circ}$C');
    plot(z_star_Cell{2,1},rad2deg(tm0_Cell{3,1}),'x',Color='black',DisplayName='$T_{\rm{o}}^*=300^{\circ}$C');
    yticks([0 10 20 30 40])
    yticklabels({'0$^\circ$' '10$^\circ$' '20$^\circ$' '30$^\circ$' '40$^\circ$'})
    ylabel('$\theta_{\mathrm{m}}$',FontSize=fsAxi);
    xlabel('$z^*$ (m)',FontSize=fsAxi); 
    xlim([0 0.1]);
    lgd = legend(Location='best',Box='off',FontSize=fsLgd)  

    nexttile
    hold on; grid on; box on
    plot(z_star_Cell{1,2},rad2deg(tm0_Cell{1,2}),'o',Color='black',DisplayName='$T_{\rm{o}}^*=60^{\circ}$C');
    plot(z_star_Cell{2,2},rad2deg(tm0_Cell{2,2}),'^',Color='black',DisplayName='$T_{\rm{o}}^*=200^{\circ}$C');
    plot(z_star_Cell{2,2},rad2deg(tm0_Cell{3,2}),'x',Color='black',DisplayName='$T_{\rm{o}}^*=300^{\circ}$C');
    yticks([0 10 20 30 40 50])
    yticklabels({'0$^\circ$' '10$^\circ$' '20$^\circ$' '30$^\circ$' '40$^\circ$' '50$^\circ$'})
    xlabel('$z^*$ (m)',FontSize=fsAxi); 
    xlim([0 0.2]);
    lgd = legend(Location='best',Box='off',FontSize=fsLgd)  

picturewidth = 17*2; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','MarkerSize'),'MarkerSize',8)
% set(findall(hfig,'-property','LineWidth'),'LineWidth',1.2)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
movegui('south')
print(hfig,hfigname,'-dpdf','-painters')
% print(hfig,hfigname,'-dpng','-painters')


%% Plotting Section (abs p* gradient ratio VERISION 1) -------------------------------------
clear all; close all; clc; format  compact; format long

hfigname = ['Figures\abs_p_star_gradient_ratio_v2'];
filenameCell{1,1} = 'Pressure_Plots\pzplot_T60_Rgs1.18e+00_La_star1.00e-01_2025-11-27'; 
filenameCell{1,2} = 'Pressure_Plots\pzplot_T60_Rgs1.20e+00_La_star1.00e-01_2025-11-27'; 
filenameCell{1,3} = 'Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star1.00e-01_2025-11-27'; 

filenameCell{2,1} = 'Pressure_Plots\pzplot_T200_Rgs1.18e+00_La_star1.00e-01_2025-11-27'; 
filenameCell{2,2} = 'Pressure_Plots\pzplot_T200_Rgs1.20e+00_La_star1.00e-01_2025-11-27'; 
filenameCell{2,3} = 'Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star1.00e-01_2025-11-27'; 

filenameCell{3,1} = 'Pressure_Plots\pzplot_T300_Rgs1.18e+00_La_star1.00e-01_2025-11-27'; 
filenameCell{3,2} = 'Pressure_Plots\pzplot_T300_Rgs1.20e+00_La_star1.00e-01_2025-11-27'; 
filenameCell{3,3} = 'Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star1.00e-01_2025-11-27'; 

hfig = figure('DefaultAxesFontSize',12); % save the figure handle in a variable
fsLgd = 12
fsAxi = 14*1.5

ls = 1.6

% DispNameCell{1} ='$T_{\rm{o}}^*=60^{\circ}$C';
% DispNameCell{2} ='$T_{\rm{o}}^*=200^{\circ}$C';
% DispNameCell{3} ='$T_{\rm{o}}^*=300^{\circ}$C';

DispNameCell{1} ='$R_{\rm{g}}^*=1.18$ mm, $\theta_{\mathrm{m}} = 0$ (case C$_0$)';
DispNameCell{2} ='$R_{\rm{g}}^*=1.20$ mm, $\theta_{\mathrm{m}} = 0$ (case C$_0$)';
DispNameCell{3} ='$R_{\rm{g}}^*=1.22$ mm, $\theta_{\mathrm{m}} = 0$ (case C$_0$)';
DispNameCell57{1} ='$R_{\rm{g}}^*=1.18$ mm, $\theta_{\mathrm{m}} = 57^{\circ}$';
DispNameCell57{2} ='$R_{\rm{g}}^*=1.20$ mm, $\theta_{\mathrm{m}} = 57^{\circ}$';
DispNameCell57{3} ='$R_{\rm{g}}^*=1.22$ mm, $\theta_{\mathrm{m}} = 57^{\circ}$';

ColorCell = {'#0072BD', '#D95319', '#77AC30'} 

% tcl = tiledlayout(3,2,"Padding","tight")
for iRg = 1:3
        load(filenameCell{1,iRg})
        ratioValue(1) = ratio_caseAC;
        ratioValue57(1) = ratio_caseAC_57;
        load(filenameCell{2,iRg})
        ratioValue(2) = ratio_caseAC;
        ratioValue57(2) = ratio_caseAC_57;
        load(filenameCell{3,iRg})
        ratioValue(3) = ratio_caseAC;
        ratioValue57(3) = ratio_caseAC_57;

    hold on; grid on; box on
    plot([60 200 300],ratioValue,'o',Color=ColorCell{iRg},LineWidth=ls,DisplayName=DispNameCell{iRg});
    plot([60 200 300],ratioValue57,'*',Color=ColorCell{iRg},LineWidth=ls,DisplayName=DispNameCell57{iRg});
    xticks([60 200 300])
    xticklabels({'60' '200' '300'})
    
    lgd = legend(Location='best',Box='off',FontSize=fsLgd)    
end 

xlabel('$T_{\rm{o}}^*$ ($^\circ$C)',FontSize=fsAxi)
title('$1/f_{\mathrm{ratio}}$',FontSize=fsAxi)

xlim([60 300]);
% ylim([0.5 2.5]);
picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','MarkerSize'),'MarkerSize',8)
% set(findall(hfig,'-property','LineWidth'),'LineWidth',1.2)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
movegui('south')
print(hfig,hfigname,'-dpdf','-painters')
print(hfig,hfigname,'-dpng','-painters')

hfigname2 = ['Figures\abs_p_star_gradient_ratio_inverse_v2'];
hfig2 = figure('DefaultAxesFontSize',12); % save the figure handle in a variable
fsLgd = 11
fsAxi = 14

ls = 1.6
% tcl = tiledlayout(3,2,"Padding","tight")
for iRg = 1:3
        load(filenameCell{1,iRg})
        ratioValue(1) = ratio_caseAC;
        ratioValue57(1) = ratio_caseAC_57;
        load(filenameCell{2,iRg})
        ratioValue(2) = ratio_caseAC;
        ratioValue57(2) = ratio_caseAC_57;
        load(filenameCell{3,iRg})
        ratioValue(3) = ratio_caseAC;
        ratioValue57(3) = ratio_caseAC_57;

    hold on; grid on; box on
    plot([60 200 300],1./ratioValue,'o',Color=ColorCell{iRg},LineWidth=ls,DisplayName=DispNameCell{iRg});
    plot([60 200 300],1./ratioValue57,'*',Color=ColorCell{iRg},LineWidth=ls,DisplayName=DispNameCell57{iRg});
    xticks([60 200 300])
    xticklabels({'60' '200' '300'})
    
    lgd = legend(Location='best',Box='off',FontSize=fsLgd)    
end 

xlabel('$T_{\rm{o}}^*$ ($^\circ$C)',FontSize=fsAxi)
ylabel('$f_{\mathrm{ratio}}$',FontSize=fsAxi)

xlim([60 300]);
% ylim([0.5 2.5]);
picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig2,'-property','MarkerSize'),'MarkerSize',8)
% set(findall(hfig,'-property','LineWidth'),'LineWidth',1.2)
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig2,'-property','Interpreter'),'Interpreter','latex')
set(findall(hfig2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
movegui('south')
print(hfig2,hfigname2,'-dpdf','-painters')
print(hfig2,hfigname2,'-dpng','-painters')


%% Velocity field
close all % dont clear all, this section will follow previous section

% for local analysis plots use:
% load('Pressure_Plots\pzplot_T60_Rgs1.18e+00_La_star1.00e-01_2025-03-11')
% load('Pressure_Plots\pzplot_T200_Rgs1.18e+00_La_star1.00e-01_2025-03-11')
% load('Pressure_Plots\pzplot_T300_Rgs1.18e+00_La_star1.00e-01_2025-03-11')

% load('Pressure_Plots\pzplot_T60_Rgs1.20e+00_La_star1.00e-01_2025-03-11')
% load('Pressure_Plots\pzplot_T200_Rgs1.20e+00_La_star1.00e-01_2025-03-11')
% load('Pressure_Plots\pzplot_T300_Rgs1.20e+00_La_star1.00e-01_2025-03-11')

% load('Pressure_Plots\pzplot_T60_Rgs1.22e+00_La_star1.00e-01_2025-03-11')
% load('Pressure_Plots\pzplot_T200_Rgs1.22e+00_La_star1.00e-01_2025-03-11')
% load('Pressure_Plots\pzplot_T300_Rgs1.22e+00_La_star1.00e-01_2025-03-11')


% ig = 3  ; % Grooved depth value
% iz = 80 ; % La* location
ig
iz
La_star(iz)
i  = 177; % \theta_m^0 location 15=0deg 177=57deg max @ 200 or 201

dxidz_v_N = Qv0N_C0(ig,iz)./tQv0N(ig,:);
dxidz_l_N = dxidz_v_N .* tXi0N(ig,:);
dxidz_v_A = Qv0A_C0(ig,iz)./tQv0A(ig,1);

meMax = 0.005 ; % max mesh element size
meMin = 0.0005; % min mesh element size

ModelU = ModelingFunction(tp,tg,RL,RgArr(ig),Rm0(i),xm(i));
% pdegplot(ModelU,EdgeLabels='on',VertexLabels='on',FaceLabels='on') % EdgeLabels for 2-D
ModelU.Mesh = generateMesh(ModelU,'Hmax',meMax,...
    'Hedge',{[6 7],meMin*0.1},'Hvertex',{[4],meMin*0.01});
% errorcheck = XiNumericalSolver(ModelU,rho,mu,tXi0N(ig,i),'rootFinding');
resultU = XiNumericalSolver(ModelU,rho,mu,tXi0N(ig,i),'velocityField');

% Domain Velocity Field % -------------------------------------------------
Nf1 = findNodes(ModelU.Mesh,"region","Face",1);
Nf2 = findNodes(ModelU.Mesh,"region","Face",2);
vz_star = zeros(length(resultU.NodalSolution),1);
vz_star(Nf1) =  resultU.NodalSolution(Nf1) .* dxidz_v_N(i)    .* V_sv_star(iz); 
vz_star(Nf2) =  resultU.NodalSolution(Nf2) .* dxidz_v_N(i)/mu .* V_sl_star(iz); 

% Mach and Reynolds number calculation 
[A,AE] = area(ModelU.Mesh); % Obtain area of each element (or get the dS):
E1_Idx = findElements(ModelU.Mesh,"region","Face",1);
E2_Idx = findElements(ModelU.Mesh,"region","Face",2);
Avapor1  = sum(AE(E1_Idx)) * P_star^2;
Aliquid1 = sum(AE(E2_Idx)) * P_star^2;
Avapor2 = tp/2*RL_star^2/2 ;
Aliquid2 = tg/2*(Rg_star(ig)^2- RL_star^2)/2;

v_zv_star_avg1 =     Qv0N_C1(ig,iz)*V_sv_star(iz)*P_star^2/Avapor1
v_zl_star_avg1 = rho*Qv0N_C1(ig,iz)*V_sv_star(iz)*P_star^2/Aliquid1
v_zv_star_avg2 = mean(vz_star(Nf1))
v_zl_star_avg2 = mean(vz_star(Nf2))

D_hydro = 4*Aliquid1./(2*(Rg_star(ig)-RL_star) + tg*(Rg_star(ig)+RL_star));
Ma = [      v_zv_star_avg1/1540          ,       v_zl_star_avg1/1540        ]
Re = [rho_v*v_zv_star_avg1*2*RL_star/mu_v, rho_l*v_zl_star_avg1*D_hydro/mu_l]

Mexy = resultU.Mesh.Nodes(:,findNodes(resultU.Mesh,"region","Edge",[7]));


% Zoom-in Domain Velocity Field% ------------------------------------------
% vapor phase:
girdMesh = 1000

[gdR_v,gdT_v] = meshgrid(linspace(0,RgArr(ig),girdMesh),linspace(0,tp/2,girdMesh));
gdX_v = gdR_v.*cos(gdT_v); gdY_v = gdR_v.*sin(gdT_v);
gdpt_v = [gdX_v(:),gdY_v(:)]';
gd_dis1  = sqrt(gdpt_v(1,:).^2+gdpt_v(2,:).^2);
gd_dis2  = sqrt((gdpt_v(1,:)-xm(i)).^2+gdpt_v(2,:).^2);
gd_ang  = atan(gdpt_v(2,:)./gdpt_v(1,:));
gd_crit_1 = (gd_dis1 > RL & gd_dis2 >= Rm0(i) & gd_ang <= tg/2);
v_zv_Zoom = interpolateSolution(resultU,gdpt_v);
v_zv_star_Zoom = v_zv_Zoom.* dxidz_v_N(i) .* V_sv_star(iz);
v_zv_star_Zoom(gd_crit_1) = NaN;
v_zv_star_Zoom = reshape(v_zv_star_Zoom,[size(gdX_v)]);

% % triple contract point phase:
% gdW_size = 0.05
% gdX_range = linspace(RL*cos(tg/2)-gdW_size,RL*cos(tg/2)+gdW_size,800); 
% gdY_range = linspace(RL*sin(tg/2)-gdW_size,RL*sin(tg/2)+gdW_size,800);
% [gdX,gdY] = meshgrid(gdX_range,gdY_range);
% gdpt = [gdX(:),gdY(:)]';
% gd_dis  = sqrt((gdpt(1,:)-xm(i)).^2+gdpt(2,:).^2);
% gd_ang  = atan(gdpt(2,:)./gdpt(1,:));
% gd_crit_v = (gd_dis <  Rm0(i) | gd_ang >  tg/2);
% gd_crit_l = (gd_dis >= Rm0(i) & gd_ang <= tg/2);

% V_z_Zoom = interpolateSolution(resultU,gdpt);
% v_z_star_Zoom = zeros(size(V_z_Zoom));
% v_z_star_Zoom(gd_crit_v) = V_z_Zoom(gd_crit_v).* dxidz_v_N(i)   .* V_sv_star(iz);
% v_z_star_Zoom(gd_crit_l) = V_z_Zoom(gd_crit_l).* dxidz_v_N(i)/mu.* V_sl_star(iz);
% v_z_star_Zoom = reshape(v_z_star_Zoom,size(gdX));

% v_zv_star_Zoom = reshape(v_zv_star_Zoom,size(gdX_v));

% liquid phase:
[gdR_l,gdT_l] = meshgrid(linspace(RL,RgArr(ig),girdMesh),linspace(0,tg/2-tg/2/girdMesh,girdMesh));
gdX_l = gdR_l.*cos(gdT_l); gdY_l = gdR_l.*sin(gdT_l);
gdpt_l = [gdX_l(:),gdY_l(:)]';
gdRm_l  = sqrt((gdpt_l(1,:)-xm(i)).^2+gdpt_l(2,:).^2);
v_zl_Zoom = interpolateSolution(resultU,gdpt_l);
v_zl_star_Zoom = v_zl_Zoom.* dxidz_v_N(i)/mu .* V_sl_star(iz);
v_zl_star_Zoom(gdRm_l < Rm0(i)) = NaN;
v_zl_star_Zoom = reshape(v_zl_star_Zoom,[size(gdX_l)]);

timestr = datetime('2025-03-12','Format','yyyy-MM-dd');
filename = sprintf('localAnalysis_T%i_Rgs%1.2i_La_star%1.2i_theta%2.0f_%s.mat',T,Rg_star(ig)*1E3,La_star(iz),rad2deg(tm0(i)),timestr)
save(['Velocity_Field_Plots\' filename]) % Saving Xi Numerical Results

%  Plotting Section (Velocity Field VERISION 1) ---------------------------------
hfig = figure('DefaultAxesFontSize',11);% hfig.Position = [0 0 1200 1200*0.75];
fsLgd = 12
fsAxi = 16

tiledlayout(1,6,"Padding","tight")
nexttile([1 5])
% pdeplot(resultU.Mesh,XYData=vz_star,Contour="on",Levels=32,ColorBar="on"); hold on
contourf(gdX_v,gdY_v,v_zv_star_Zoom,32); hold on
box on; axis equal; % plot the solution
xlim([0 4.3]);ylim([0,RL*sin(tp/2)]);
set(gca,'YTick',[])
colormap(gca,'parula');
xlabel('$r^*$',Fontsize=fsAxi,Interpreter='latex')
cb = colorbar(Location='northoutside');
ylabel(cb,'$v^*_{\rm{v}}$~[m/s]',FontSize=fsAxi,Rotation=0,Interpreter='latex')

nexttile
contourf(gdX_l,gdY_l,v_zl_star_Zoom,16); hold on
box on; axis equal; % plot the solution
xlabel('$r^*$',Fontsize=fsAxi,Interpreter='latex')
colormap('turbo');
xlim([4.1 4.5])
set(gca,'YTick',[])
cb = colorbar(Location='northoutside');
ylabel(cb,'$v^*_{\rm{l}}$~[m/s]',FontSize=fsAxi,Rotation=0,Interpreter='latex')
% title('Liquid Phase')

% Zoom-in Domain Velocity Field% 
% subplot(2,2,3);
% contourf(gdX,gdY,v_z_star_Zoom,16); hold on
% % plot(Mexy(1,:),Mexy(2,:),'--',Linewidth=0.8,Color='black') % plot meniscus
% % xlim([min(gdX,[],"all"),max(gdX,[],"all")])
% % ylim([min(gdY,[],"all"),max(gdY,[],"all")])
% axis equal;
% colormap('parula');colorbar;
% title('Triple Contract Point')

% title([ 'T = ',num2str(T),'^\circC;    ',...
%         'Rg* = ',num2str(Rg_star(ig)),'m;    ',...
%         'La* = ',num2str(La_star(iz)),'m;    ',...
%         '\theta_m^0 = ',num2str(rad2deg(tm0(i))),'^\circ;    ',...
%         'Ma_v= ',num2str(Ma(1)),', Re_v = ',num2str(Re(1)),';    ',...])
%         'Ma_l= ',num2str(Ma(2)),', Re_l = ',num2str(Re(2))])

%--------------------------------------------------------------------------

% figure % 3D plot
% pdeplot(resultU.Mesh,XYData=vz_star,ZData=vz_star,ColorBar="off");
% title('v_z^0*'); box on; grid on;  % plot the solution
% colormap(gca,'hsv');colorbar('Location', 'southoutside');
% figure % 3D plot
% pdeplot(resultU.Mesh,XYData=resultU.NodalSolution,ZData=resultU.NodalSolution,ColorBar="off");
% title('V_z^0'); box on; grid on;  % plot the solution
% colormap(gca,'hsv');colorbar('Location', 'southoutside');

%% Plotting Section (Velocity Field VERISION 2) ---------------------------------
clear all; close all; clc; format  compact; format long g

% load('Velocity_Field_Plots\localAnalysis_T60_Rgs1.22e+00_La_star1.00e-01_theta 0_2025-03-12')
load('Velocity_Field_Plots\localAnalysis_T60_Rgs1.22e+00_La_star1.00e-01_theta57_2025-03-12')
% load('Velocity_Field_Plots\localAnalysis_T200_Rgs1.22e+00_La_star1.00e-01_theta 0_2025-03-12')
% load('Velocity_Field_Plots\localAnalysis_T200_Rgs1.22e+00_La_star1.00e-01_theta57_2025-03-12')
% load('Velocity_Field_Plots\localAnalysis_T300_Rgs1.22e+00_La_star1.00e-01_theta 0_2025-03-12')
% load('Velocity_Field_Plots\localAnalysis_T300_Rgs1.22e+00_La_star1.00e-01_theta57_2025-03-12')

% load('Velocity_Field_Plots\localAnalysis_T60_Rgs1.18e+00_La_star1.00e-01_theta 0_2025-03-12')
% load('Velocity_Field_Plots\localAnalysis_T60_Rgs1.18e+00_La_star1.00e-01_theta57_2025-03-12')

% filenameCell{1} = 'Velocity_Field_Plots\localAnalysis_T60_Rgs1.22e+00_La_star1.00e-01_theta 0_2025-03-12'; 
% filenameCell{2} = 'Velocity_Field_Plots\localAnalysis_T60_Rgs1.22e+00_La_star1.00e-01_theta57_2025-03-12'; 
% filenameCell{3} = 'Velocity_Field_Plots\localAnalysis_T200_Rgs1.22e+00_La_star1.00e-01_theta 0_2025-03-12'; 
% filenameCell{4} = 'Velocity_Field_Plots\localAnalysis_T200_Rgs1.22e+00_La_star1.00e-01_theta57_2025-03-12'; 
% filenameCell{5} = 'Velocity_Field_Plots\localAnalysis_T300_Rgs1.22e+00_La_star1.00e-01_theta 0_2025-03-12'; 
% filenameCell{6} = 'Velocity_Field_Plots\localAnalysis_T300_Rgs1.22e+00_La_star1.00e-01_theta57_2025-03-12'; 

hfig3 = figure('DefaultAxesFontSize',11);% hfig.Position = [0 0 1200 1200*0.75];
fsLgd = 16
fsAxi = 16
lwZro = 3

% tiledlayout(6,5,"Padding","tight")
tiledlayout(1,5,"Padding","tight")
% for i_vf = 1:6
% load(filenameCell{i_vf})
lgdText =  sprintf('$T_{\\rm{o}}^* = %i^{\\circ}\\rm{C}, \\theta_{\\rm{m}}^0 = %2.0f^{\\circ}$',T,rad2deg(tm0(i)))

nexttile([1 4])
% pdeplot(resultU.Mesh,XYData=vz_star,Contour="on",Levels=32,ColorBar="on"); hold on
contourf(gdX_v,gdY_v,v_zv_star_Zoom,32,HandleVisibility="off"); hold on
contour(gdX_v,gdY_v,v_zv_star_Zoom, [0, 0],':w', 'LineWidth',lwZro,HandleVisibility="off");
box on; axis equal; % plot the solution
xlim([0 4.3]);ylim([0,RL*sin(tp/2)]);
set(gca,'YTick',[])
colormap(gca,'parula');
xlabel('$r$',Fontsize=fsAxi,Interpreter='latex')
cb = colorbar(Location='northoutside');
ylabel(cb,'$v^*_{\rm{v}}$~[m/s]',FontSize=fsAxi,Rotation=0,Interpreter='latex')
lgd = legend(Location="northwest",Box='off')
title(lgd,lgdText,FontSize=fsLgd)

nexttile([1 1])
contourf(gdX_l,gdY_l,v_zl_star_Zoom,12); hold on
contour(gdX_l,gdY_l,v_zl_star_Zoom, [0, 0],':w', 'LineWidth',lwZro);
box on; axis equal; % plot the solution
xlabel('$r$',Fontsize=fsAxi,Interpreter='latex')
colormap('turbo');
xlim([4.1 4.7])
set(gca,'YTick',[])
cb = colorbar(Location='northoutside');
if T == 60
cb.Ruler.Exponent=0;                   % set the desired exponent
% cb.Ruler.TickLabelFormat='%0.3f';       % fix up ugly default %g formatting
end 
ylabel(cb,'$v^*_{\rm{l}}$~[m/s]',FontSize=fsAxi,Rotation=0,Interpreter='latex')
% end 

hfigname = ['Figures\vz_star_Field' sprintf('_T%i',T) '_Rgstr' num2str(Rg_star(ig)*100000) '_Lastr' num2str(La_star(iz)*100) sprintf('_theta%2.0f',rad2deg(tm0(i)))];
% hfigname = ['Figures\vz_star_Field'];
picturewidth = 17*2; % set this parameter and keep it forever
hw_ratio = 0.2; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig3,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig3,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig3,'Position');
set(hfig3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig3,hfigname,'-dpng','-painters')



%% Local Analysis ---------------------------------------------------------
% Velocity along Meniscus-------------------------------------------------- plot--
% pdeplot(resultU.Mesh,XYData=resultU.NodalSolution)
EdgeID = findNodes(resultU.Mesh,"region","Edge",[6 7]); % nodes id of the meniscus
Edgexy = resultU.Mesh.Nodes(:,EdgeID)    ; % x and y values of the Edge 6
thetaN  = atan(Edgexy(2,:)./Edgexy(1,:))   ;  % Edgexy6(2,:) is the y valuse of nodes on Edge 6 (Wall)
% Nondimentional V^0 and V^00 plot
Uz0  = interpolateSolution(resultU,Edgexy)';
% Andrew's solution section
[thetaA,Vzv,vzv] = A_LocalAnalysisPlots(mu,tp,tg,RL,RgArr(ig),P_star,...
                    h_lv,rho_v,rho_l,tXi0A(ig,1),tXi0A(ig,2),N_An,Ng,tm0(i),...
                    0,0,0,0,0,'C0');

% Meniscus Velocity Plot
% v_zv_starN_menis =   Uz0 * dxidz_v_N(i) * V_sv_star(iz);
% v_zv_starA_menis =   dxidz_v_A * Vzv * V_sv_star(iz);
% % v_zv_starN_menis2 = -vz_star(EdgeID)'*mu./tXi0N(ig,i)*V_sv_star(iz)/V_sl_star(iz);

% fs = 16;
% figure; 
% plot(thetaN,v_zv_starN_menis,'.',MarkerSize=10,DisplayName='numerical code');hold on
% plot(thetaA,v_zv_starA_menis,'-',LineWidth=1.75,DisplayName='Andrews case 00')
% % plot(thetaN,v_zv_starN_menis2,'o',LineWidth=1.75,DisplayName='numerical code2')
% box on; axis square
% xlabel('$\theta$',Interpreter='latex',FontSize=fs);
% ylabel('$v_{z,v}^*$ [m/s]',Interpreter='latex',FontSize=fs)
% xlim([0,tp/2])
% legend(Location="southeast")
% title([ 'T = ',num2str(T),'^\circC;  ',...
%         'Rg* = ',num2str(Rg_star(ig)),'m;  '],...
%         ['La* = ',num2str(La_star(iz)),'m;  ',...
%         '\theta_m^0 = ',num2str(rad2deg(tm0(i))),'^\circ; '])


fs = 16
ms = 10

% Below fix the issue that evaluateGradient using results in Liquid phase 
EdegLogical = thetaN < tg/2; 
Edgexy_Corrected(1,:) = Edgexy(1,:)-meMin.*0.45.*cos(thetaN).*EdegLogical;
Edgexy_Corrected(2,:) = Edgexy(2,:)-meMin.*0.45.*cos(thetaN).*EdegLogical;

[Ugradx,Ugrady] = evaluateGradient(resultU,Edgexy_Corrected);
Ugradr = Ugradx'.* cos(thetaN) + Ugrady'.* sin(thetaN);

% Andrew's solution section
[thetaA,Vzv,dVzvdr] = A_LocalAnalysisPlots(mu,tp,tg,RL,RgArr(ig),P_star,...
                    h_lv,rho_v,rho_l,tXi0A(ig,1),tXi0A(ig,2),N_An,Ng,tm0(i),...
                    0,0,0,0,0,'C0');

if tm0(i) < pi/3
    Lamig = 0.9;
else 
    Lamig = 0.8;
end 

funLam = @(lambda) tan(lambda*(pi+tm0(i))) - mu*tan(lambda*(tm0(i)-pi/2)); % function
lambda = fzero(funLam,Lamig)
localAnasMult = 0.75*abs(Uz0(end-1))/((tg/2-thetaN(end-1)).^lambda)

% Local analysis Plot section
close all

timestr = datetime('2025-03-11','Format','yyyy-MM-dd');
filename = sprintf('localAnalysis_T%i_Rgs%1.2i_La_star%1.2i_%s.mat',T,Rg_star(ig)*1E3,La_star(iz),timestr)
save(['Local_Analysis_Plot\' filename]) % Saving Xi Numerical Results

hfig = figure('DefaultAxesFontSize',12); % save the figure handle in a variable
fsLgd = 14
fsAxi = 14

tiledlayout(1,3,"Padding","tight")

% subplot(1,3,1)% -----------------------------------------------------------
nexttile
loglog(tg/2-thetaN,abs(Uz0),'.',MarkerSize=ms,...
    DisplayName='Numerical');hold on
if tm0(i) < 0.5
loglog(tg/2-thetaA,abs(Vzv),'-',Color='red',...
    DisplayName='Analytical');
end 
loglog(tg/2-thetaA,localAnasMult*(tg/2-thetaA).^lambda,'black',...
    DisplayName='Slope of ($\theta_g/2-\theta)^\lambda$')
box on; axis square
xlim([10^-6 10^-3])
xlabel('$\theta_g/2 - \theta$[rad]',FontSize=fsAxi)
ylabel('$|V_{z,v}|$',FontSize=fsAxi)
lgd = legend(Location='northwest',FontSize=fsLgd,Box='off')
% legend(Location="southeast",FontSize=fsLgd,Box='off')
title(lgd,['$\lambda = ',num2str(lambda) '$'])


% subplot(1,3,2)% -----------------------------------------------------------
nexttile
% ArcLenA = (tp/2-thetaA)*RL;
% ArcLenFunc = @(targ) (targ <= (tp/2 - tg/2)).*targ.*RL + ...
%                      (targ > (tp/2 - tg/2)).*(((tp/2-tg/2)*RL)+...
%                      (asin(xm(i)./Rm0(i).*sin(tp/2-targ))+tp/2-targ).*Rm0(i));
% ArcLenN = ArcLenFunc(tp/2-thetaN);
% 
% plot(tp/2-thetaN,ArcLenN); hold on
% plot(tp/2-thetaA,ArcLenA)

plot(thetaN, Uz0,'.',MarkerSize=ms,DisplayName='Numerical');hold on
plot(thetaA,-Vzv,'-',LineWidth=1.75,DisplayName='Analytical')
% plot(thetaN,v_zv_starN_menis2,'o',LineWidth=1.75,DisplayName='numerical code2')
box on; axis square
xlabel('$\theta$',FontSize=fsAxi);
ylabel('$V_{z,v}^0$',FontSize=fsAxi)
xlim([0,0.12])
% legend(Location="northwest",FontSize=fsLgd,Box='off')
legend(Location="northeast",FontSize=fsLgd,Box='off')
% legend(Location="southwest",FontSize=fsLgd,Box='off')


% subplot(1,3,3)% -----------------------------------------------------------
nexttile
plot(thetaN, Ugradr,'.',MarkerSize=ms,DisplayName='Numerical');hold on
plot(thetaA,-dVzvdr,'-',LineWidth=1.75,DisplayName='Analytical')
% plot(thetaN,Ugradr./dVzvdr,'o',LineWidth=1.75,DisplayName='numerical code2')
box on; axis square
xlabel('$\theta$',FontSize=fsAxi);
ylabel('$\partial V_{z,v}^0/\partial r$',FontSize=fsAxi)
xlim([0,0.12])
legend(Location="southwest",FontSize=fsLgd,Box='off')

hfigname = ['Figures\LocalAnalysis' sprintf('_T%i',T) '_Rgstr' num2str(Rg_star(ig)*100000) '_Lastr' num2str(La_star(iz)*100) sprintf('_theta%2.0f',rad2deg(tm0(i)))];

picturewidth = 17*2; % set this parameter and keep it forever
hw_ratio = 0.40; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
% print(hfig,hfigname,'-dpng','-painters')

%% Local Analysis Verision II---------------------------------------------------------
clear all; close all; clc; format  compact; format long

filenameCell{1} = 'Local_Analysis_Plot\localAnalysis_T60_Rgs1.22e+00_La_star1.00e-01_2025-03-11'; 
filenameCell{2} = 'Local_Analysis_Plot\localAnalysis_T200_Rgs1.22e+00_La_star1.00e-01_2025-03-11'; 
filenameCell{3} = 'Local_Analysis_Plot\localAnalysis_T300_Rgs1.22e+00_La_star1.00e-01_2025-03-11'; 

% filenameCell{1} = 'Local_Analysis_Plot\localAnalysis_T60_Rgs1.20e+00_La_star1.00e-01_2025-03-11'; 
% filenameCell{2} = 'Local_Analysis_Plot\localAnalysis_T200_Rgs1.20e+00_La_star1.00e-01_2025-03-11'; 
% filenameCell{3} = 'Local_Analysis_Plot\localAnalysis_T300_Rgs1.20e+00_La_star1.00e-01_2025-03-11'; 

% filenameCell{1} = 'Local_Analysis_Plot\localAnalysis_T60_Rgs1.18e+00_La_star1.00e-01_2025-03-11'; 
% filenameCell{2} = 'Local_Analysis_Plot\localAnalysis_T200_Rgs1.18e+00_La_star1.00e-01_2025-03-11'; 
% filenameCell{3} = 'Local_Analysis_Plot\localAnalysis_T300_Rgs1.18e+00_La_star1.00e-01_2025-03-11'; 

hfig2 = figure('DefaultAxesFontSize',10); % save the figure handle in a variable
fsLgd = 12
fsAxi = 14

tiledlayout(2,3,"Padding","tight")

nexttile
    load(filenameCell{1});
    plot(thetaN, Uz0,'.',MarkerSize=ms,DisplayName='Numerical');hold on
    plot(thetaA,-Vzv,'-',LineWidth=1.75,DisplayName='Analytical',Color='#D95319')
    % plot(thetaN,v_zv_starN_menis2,'o',LineWidth=1.75,DisplayName='numerical code2')
    box on; %axis square
    xlabel('$\theta$ (rad)',FontSize=fsAxi);
    ylabel('$V_{z,\rm{v}}^0$',FontSize=fsAxi)
    xlim([0,0.12])
    lgd = legend(Box='off',Location='west',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=60^{\circ}$C')

nexttile
    load(filenameCell{2});
    plot(thetaN, Uz0,'.',MarkerSize=ms,DisplayName='Numerical');hold on
    plot(thetaA,-Vzv,'-',LineWidth=1.75,DisplayName='Analytical',Color='#D95319')
    % plot(thetaN,v_zv_starN_menis2,'o',LineWidth=1.75,DisplayName='numerical code2')
    box on;% axis square
    xlabel('$\theta$ (rad)',FontSize=fsAxi);
    % ylabel('$V_{z,v}^0$',FontSize=fsAxi)
    xlim([0,0.12])
    ylim([-0.25 0])
    lgd = legend(Box='off',Location='east',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=200^{\circ}$C')

nexttile
    load(filenameCell{3});
    plot(thetaN, Uz0,'.',MarkerSize=ms,DisplayName='Numerical');hold on
    plot(thetaA,-Vzv,'-',LineWidth=1.75,DisplayName='Analytical',Color='#D95319')
    % plot(thetaN,v_zv_starN_menis2,'o',LineWidth=1.75,DisplayName='numerical code2')
    box on;% axis square
    xlabel('$\theta$ (rad)',FontSize=fsAxi);
    % ylabel('$V_{z,v}^0$',FontSize=fsAxi)
    xlim([0,0.12])
    lgd = legend(Box='off',Location='east',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=300^{\circ}$C')

% subplot(1,3,1)% -----------------------------------------------------------
nexttile
    load(filenameCell{1});
    plot(thetaN, Ugradr,'.',MarkerSize=ms,DisplayName='Numerical');hold on
    plot(thetaA,-dVzvdr,'-',LineWidth=1.75,DisplayName='Analytical',Color='#D95319')
    % plot(thetaN,Ugradr./dVzvdr,'o',LineWidth=1.75,DisplayName='numerical code2')
    box on;% axis square
    xlabel('$\theta$ (rad)',FontSize=fsAxi);
    ylabel('$\partial V_{z,\rm{v}}^0/\partial r$',FontSize=fsAxi)
    xlim([0,0.12])
    lgd = legend(Box='off',Location='west',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=60^{\circ}$C')
    
    % subplot(1,3,2)% -----------------------------------------------------------
nexttile
    load(filenameCell{2})
    plot(thetaN, Ugradr,'.',MarkerSize=ms,DisplayName='Numerical');hold on
    plot(thetaA,-dVzvdr,'-',LineWidth=1.75,DisplayName='Analytical',Color='#D95319')
    % plot(thetaN,Ugradr./dVzvdr,'o',LineWidth=1.75,DisplayName='numerical code2')
    box on; % axis square
    xlabel('$\theta$ (rad)',FontSize=fsAxi);
    % ylabel('$\partial V_{z,v}^0/\partial r$',FontSize=fsAxi)
    xlim([0,0.12])
    lgd = legend(Box='off',Location='west',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=200^{\circ}$C')

    % subplot(1,3,3)% -----------------------------------------------------------
nexttile
    load(filenameCell{3})
    plot(thetaN, Ugradr,'.',MarkerSize=ms,DisplayName='Numerical');hold on
    plot(thetaA,-dVzvdr,'-',LineWidth=1.75,DisplayName='Analytical',Color='#D95319')
    % plot(thetaN,Ugradr./dVzvdr,'o',LineWidth=1.75,DisplayName='numerical code2')
    box on; % axis square
    xlabel('$\theta$ (rad)',FontSize=fsAxi);
    % ylabel('$\partial V_{z,v}^0/\partial r$',FontSize=fsAxi)
    xlim([0,0.12])
    lgd = legend(Box='off',Location='west',FontSize=fsLgd);
    title(lgd,'$T_{\rm{o}}^*=300^{\circ}$C')

hfigname = ['Local_Analysis_Plot\LocalAnalysis_' 'allT' '_Rgstr' num2str(Rg_star(ig)*100000) '_Lastr' num2str(La_star(iz)*100)];

picturewidth = 17*1.75; % set this parameter and keep it forever
hw_ratio = 0.6; % feel free to play with this ratio
% set(findall(hfig2,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig2,'Units','centimeters','Position',[1 1 picturewidth hw_ratio*picturewidth])
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig2,hfigname,'-dpdf','-painters')
print(hfig2,hfigname,'-dpng','-painters')

%% Effects of Gravity
clear all; close all; clc; format  compact; format long
load('XiRgs1.18e+00-1.22e+00_T60_2025-01-07.mat')
% load('XiRgs1.18e+00-1.22e+00_T200_2025-01-09.mat')
% load('XiRgs1.18e+00-1.22e+00_T300_2025-01-07.mat')

clear Model meMax meMin cl  earance filename timestr NIST
clear tm0_max tm0_min Rg epsilon V_sv_star V_sl_star
clear tQv00 tQv01 tQl00 tQl01 tXi00 tXi01 tXi00_ig tXi0N_ig

% Solver Settings ---------------------------------------------------------
N_pt = 200;  % numbers of interpulations points

% Remove the NaN term
  tm0 = tm0(1:200);
  Rm0 = Rm0(1:200);
tXi0N = tXi0N(:,1:200);
tQv0N = tQv0N(:,1:200);

% Cross-Sectional Variables -----------------------------------------------
Dh = 2*tg*(Rg_star.^2 - RL_star^2)./...
    (2*(Rg_star - RL_star) + tg.*(Rg_star+RL_star) );

phi_l = Ng*tg/(2*pi);

% Wick permeability per Sparrow et al. (1963)
SparrowRatio = [ 0.7   0.8   0.85  0.9   0.95]; 
SparrowValue = [18.08 16.42 15.32 14.35 15.06;...
                15.56 14.42 14.26 14.98 17.62;...
                14.52 14.33 14.91 16.36 19.17;...
                14.26 14.84 15.82 17.51 20.14;...
                14.71 16.15 17.39 19.07 21.26;...
                15.52 17.28 18.52 20.05 21.88;...
                16.33 18.16 19.34 20.71 22.28;...
                17.04 18.86 19.95 21.18 22.56];
SparrowAngel = deg2rad([5 10 15 20 30 40 50 60]);
Sparrow_fRe = interp2(SparrowRatio,SparrowAngel,SparrowValue,RL_star./Rg_star,tg)
clear SparrowRatio SparrowValue SparrowAngel

K = Dh.^2 * phi_l./(2*Sparrow_fRe)
Aw_star = pi*(Rg_star.^2 - RL_star^2);


% Longitudinal Variables --------------------------------------------------
g = 9.81;
% phi = deg2rad(0.01); % tilting angle in rad, + for 
% phi =deg2rad([0 0.01 0.02 0.03]); % tilting angle in rad, + for 
% phi = linspace(-0.2673,0.525,150); % tilting angle in rad, + for 
phi = linspace(-0.2673,pi + 0.2673,300); % tilting angle in rad, + for 
phi = [phi(phi<0) 0 phi(phi<pi/2& phi>0) pi/2 phi(phi>pi/2 )];


    % \tilde \xi for Case A & B -------------------------------------------
    tXi00_CaseA = - pi*rho*RL_star^4./(8*K.*Aw_star*mu);
    
    for i_geom = 1:length(RgArr) 
        [K1_star,K2_star,K4_star,Ca_star] = A_FullMassFlowRate_CaseAB(mu,tp,tg,RL_star,RL,Rg_star(i_geom),RgArr(i_geom),P_star,tXi0A(i_geom,1),tXi0A(i_geom,2),N_An);
        tXi00_CaseB(i_geom) = (((pi*rho*RL_star^4)/(16*Ng*mu))-K4_star*RL_star/2) / (K1_star+K2_star-Ca_star*K4_star);
        clear K1_star K2_star K4_star Ca_star
    end
    
    % ---------------------------------------------------------------------

La_star = 0.100; % This has changed to 5x larger
Le_star = 0.100;
Lc_star = 0.100;
% La_star = linspace(0.1,1,100);
Lt_star = 0.1 + 0.1 + 0.1; % Total pipe length
% La_star = linspace(Lt_star*0.975,Lt_star*0.025,100)

Leff_star = Le_star/2 + La_star + Lc_star/2;

epsilon =  P_star./La_star     ; % small parameter
Dp_cmax_star = sigma/Rp_star;
   Ca_hat = Dp_cmax_star*P_star/sigma;
V_sv_star = Dp_cmax_star*P_star*epsilon/mu_v;
V_sl_star = Dp_cmax_star*P_star*epsilon/mu_l;

% Non-dimensionalization --------------------------------------------------
  Aw = Aw_star.*K/(P_star^4);
  Le = Le_star./La_star     ;
  Lc = Lc_star./La_star     ;
Bo_v = rho_v*g.*La_star*sin(phi)/Dp_cmax_star;
Bo_l = rho_l*g.*La_star*sin(phi)/Dp_cmax_star;
d_Bo = (rho_l-rho_v)*g.*La_star*sin(phi)/Dp_cmax_star; % ΔBo

% Capillary limit calculation and plots
clc
Qv0A_C0_g = zeros(length(phi),1);
Qv0N_C0_g = zeros(length(phi),1);
Qv0N_C1_g = zeros(length(phi),1);

qc_CaseA_g = zeros(length(phi),1);
qc_CaseB_g = zeros(length(phi),1);

switch T
    case  60; tQv0_igA =  0; %0.005;%0.0009;%[0.01 1]; % 60C-1
    case 200; tQv0_igA =  0.1; %0.005;%0.0009;%[0.01 1]; % 200C-0.1
    case 300; tQv0_igA =  0.006; %0.005;%0.0009;%[0.01 1]; % 300C-0.01
end

 ig = 3 % loop counter for different domain shape

    SC.rho = rho;
    SC.mu = mu;
    SC.Ng = Ng;
    SC.RL = RL;
    SC.tg = tg;
    SC.Aw = Aw(ig);
    SC.Ca_hat = Ca_hat;

    tXi00 = tXi0A(ig,1);
    tXi01 = tXi0A(ig,2);
    tQv00 = tQv0A(ig,1);
    tQv01 = tQv0A(ig,3);


    for j = 1:length(phi)  % loop counter for different titling angle

        % Longitudinal Varaible Structure:
        SL.Le = Le;
        SL.Lc = Lc;
        SL.d_Bo = d_Bo(j);

        % Qv^0 Initial Guess ----------------------------------------------
        if j == 1
            tQv0_igN =  tQv0_igA ;
            tQv0_igN1 = tQv0_igA ;
        else
            tQv0_igA = Qv0A_C0_g(j-1); % tilde xi 0 initial guess 1
            tQv0_igN = Qv0A_C0_g(j-1);
            tQv0_igN1 = Qv0N_C1_g(j-1);
        end

        Qv0A_C0_g(j) = fzero(@(Qv0Arg) QvRootAndrew(SC,SL,tXi00,tXi01,tQv00,tQv01,Qv0Arg,'C0'),tQv0_igA);
        % Qv0N_C0_g(j) = fzero(@(Qv0Arg) ...
        %     QvRootNumerical_g(SC,SL,tm0,Rm0,tXi0N(ig,:),tQv0N(ig,:),Qv0Arg,N_pt,'C0'),tQv0_igN);
        % Qv0N_C1_g(j) = fzero(@(Qv0Arg) ...
            % QvRootNumerical_g(SC,SL,tm0,Rm0,tXi0N(ig,:),tQv0N(ig,:),Qv0Arg,N_pt,'C1'),tQv0_igN1);
    end

    clear SC SL

    % Case A & B
    qc_CaseA_g = (Dp_cmax_star - (rho_v - rho_l)*g*sin(phi)*Lt_star)./...
    (8*mu_v*Leff_star/(pi*rho_v*RL_star^4*h_lv)+...
    mu_l*Leff_star./(K(ig)*rho_l*Aw_star(ig)*h_lv));

    [K1_star,K2_star,K4_star,Ca_star] = A_FullMassFlowRate_CaseAB(mu,tp,tg,RL_star,RL,Rg_star(ig),RgArr(ig),P_star,tXi0A(ig,1),tXi0A(ig,2),N_An);
    qc_CaseB_g = (Dp_cmax_star - (rho_v - rho_l)*g*sin(phi)*Lt_star)./...
    (Leff_star*( 8*mu_v/(pi*rho_v*RL_star^4*h_lv) + (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv) )/(K1_star + K2_star - Ca_star*K4_star)*mu_l ));



tQv0N(:,1)

qc0A_C0_g = Qv0A_C0_g.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;
qc0N_C0_g = Qv0N_C0_g.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;
qc0N_C1_g = Qv0N_C1_g.*V_sv_star*P_star^2*rho_v*h_lv*2*Ng;

timestr = datetime('2025-07-17','Format','yyyy-MM-dd');
filename = sprintf('EffGCapLimRgs%1.2i-%1.2i_T%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,timestr)
save(['Effects_Gravity\' filename])  % Saving Xi Numerical Results


ratioCap0bypiby2 = qc0A_C0_g(152)/qc0A_C0_g(23) * 100

% Plotting Section (Capillary Limit VERISION 1) -------------------------------------

CArr = ["#0072BD","#D95319","#77AC30","black"];

close all; 
hfig = figure('DefaultAxesFontSize',11);
tiledlayout(1,1,"Padding","tight")

lw = 1.4

nexttile
hold on; box on
    % nameStr = '$q^{*}_{\mathrm{c,max}}$';
    % plot(phi,qc0N_C1_g,'.'  ,LineWidth = lw,Color=CArr(1),DisplayName=nameStr)
    nameStr = 'C$_0$';
    plot(rad2deg(phi),qc0A_C0_g,'-' ,LineWidth = lw,Color=CArr(2),DisplayName=nameStr)
    nameStr = 'A';
    plot(rad2deg(phi),qc_CaseA_g,'--' ,LineWidth = lw,Color=CArr(3),DisplayName=nameStr)
    nameStr = 'B';
    plot(rad2deg(phi),qc_CaseB_g,'-.' ,LineWidth = lw,Color=CArr(4),DisplayName=nameStr)
    xlim(rad2deg([-pi/6  pi*7/6]))
    xticks(rad2deg([-pi/6 0 pi/6 pi/3 pi/2 pi*2/3 pi*5/6 pi pi*7/6]))
    xticklabels({'-30°' '0°' '30°' '60°' '90°' '120°' '150°' '180°' '210°'})
    xlabel('$\phi$') % (Lt* = ',num2str(Lt_star),' m)'])
    ylabel('$q^*_{\mathrm{cap}}$ (W)')
    % title(['Capillary limit vs Adiabatic sec length at T = ',num2str(T),'°C'])
    legend(Location="south",Box="off",Orientation="horizontal")

% nexttile
%     % Angle plots!
%     p_ve_g = 1 - 8*Ng*Qv0A_C0_g'*Le./(pi*RL^4) - 1*Ng*rho*Qv0A_C0_g'*Le./(mu*Aw(ig))+d_Bo*Le;
%     delta_pc_g  = 8*Ng*Qv0A_C0_g'*Lc./(pi*RL^4) + 1*Ng*rho*Qv0A_C0_g'*Lc./(mu*Aw(ig))-d_Bo*Lc;
%     Rm0_e_g  = 1./(Ca_hat*p_ve_g);
%     Rm0_c_g  = 1./(Ca_hat*delta_pc_g);
%     hold on; box on
%     nameStr = 'condensor inlet (C$_0$)';
%     plot(phi,Rm0_c_g,'-'  ,LineWidth = lw,Color='k',DisplayName=nameStr)
%     nameStr = 'evaporator outlet (C$_0$)';
%     plot(phi,Rm0_e_g,'--' ,LineWidth = lw,Color='k',DisplayName=nameStr)
%     xlabel('$\phi$ (rad)') % (Lt* = ',num2str(Lt_star),' m)'])
%     ylabel('$R^0_{\mathrm{m}}$')
%     legend(Location="best",Box="off",Orientation="vertical")

hfigname = ['Figures\GravityEff_baseline_noC1_v1' sprintf('_T%i',T)];
% hfigname = ['Figures\GravityEff_baseline_noC1_v1' sprintf('_T%i',T)];
picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 0.4; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document
% set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
% print(hfig,hfigname,'-dpng','-painters')



%% Effects of Gravity: p* vs z* Deep Dive Study 
clear all; close all; clc; format  compact; format long

load('Effects_Gravity\EffGCapLimRgs1.18e+00-1.22e+00_T60_2025-07-17'); 
N_pt = 1000;

ig  % Grooved depth value
ip = 151
La_star
i  = 15; % \theta_m^0 location

dxidz_v_A = Qv0A_C0_g(ip)./tQv0A(ig,1);
dxidz_l_A = dxidz_v_A .* tXi0A(ig,1);

ze_star = linspace(-Le_star,              0, 50);
zc_star = linspace( La_star,La_star+Lc_star, 50);
za      = linspace(        0,             1, 50);
za_star = za.*La_star;
z_star = [ze_star za_star zc_star];


for i=i % Analytical Case C0
p_veA = 1 - 8*Ng*Qv0A_C0_g(ip)*Le/(pi*RL^4) ...
         - 1*Ng*rho*Qv0A_C0_g(ip)*Le/(mu*Aw(ig))+d_Bo(ip)*Le;
d_pcA = 8*Ng*Qv0A_C0_g(ip)*Lc/(pi*RL^4) ...
         + 1*Ng*rho*Qv0A_C0_g(ip)*Lc/(mu*Aw(ig))-d_Bo(ip)*Lc;

d_pe_starA =  p_veA .* Dp_cmax_star;
d_pc_starA =  d_pcA .* Dp_cmax_star;


% Evaporator-----------
p_ve_starA = p_sat - 8*mu_v*qc0A_C0_g(ip)./(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star)+ze_star+Le_star/2) ...
    - (rho_v*g*sin(phi(ip)))*(Le_star+ze_star);
p_le_starA = p_sat - Dp_cmax_star + mu_l*qc0A_C0_g(ip)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star)+ze_star+Le_star/2)...
    - rho_l*g*sin(phi(ip))*(Le_star+ze_star);
diffcheck_evap = p_ve_starA(end) - p_le_starA(end) - d_pe_starA 

% Adiabatic  -----------
p_vaA = (-Bo_v(ip) + dxidz_v_A) .* za + d_pe_starA/Dp_cmax_star; % pva
p_laA = (-Bo_l(ip) + dxidz_l_A) .* za;
p_va_starA = p_vaA.*Dp_cmax_star + p_le_starA(end);
p_la_starA = p_laA.*Dp_cmax_star + p_le_starA(end);
diffcheck_adia = p_va_starA(end) - p_la_starA(end) - d_pc_starA 

% Condnesor  -----------
p_vc_starA = p_va_starA(end) - 8*mu_v*qc0A_C0_g(ip)./(pi*rho_v*RL_star^4*h_lv)*...
      ( (La_star+Lc_star).*zc_star - zc_star.^2/2 - La_star.^2/2 - La_star*Lc_star)/Lc_star... 
      - rho_v*g*sin(phi(ip))*(zc_star-La_star);
p_lc_starA = p_la_starA(end) + mu_l*qc0A_C0_g(ip)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
      ( (La_star+Lc_star).*zc_star - zc_star.^2/2 - La_star.^2/2 - La_star*Lc_star)/Lc_star... 
      - rho_l*g*sin(phi(ip))*(zc_star-La_star);

p_v_starA = [p_ve_starA p_va_starA p_vc_starA];
p_l_starA = [p_le_starA p_la_starA p_lc_starA];
end 

for i=i  % Analytical  Case A
delta_p_vt_A = 8*mu_v*qc_CaseA_g(ip)*Leff_star/(pi*rho_v*RL_star^4*h_lv)+...
    rho_v*g*sin(phi(ip))*Lt_star;

p_ve_star_A = p_sat - 8*mu_v*qc_CaseA_g(ip)/(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star)+ze_star+Le_star/2) - rho_v*g*sin(phi(ip))*...
    (Le_star+ze_star);
p_va_star_A = p_sat - 8*mu_v*qc_CaseA_g(ip)/(pi*rho_v*RL_star^4*h_lv)*...
    (za_star+Le_star/2) - rho_v*g*sin(phi(ip))*(Le_star+za_star);
p_vc_star_A = p_sat - 8*mu_v*qc_CaseA_g(ip)/(pi*rho_v*RL_star^4*h_lv)*...
    (Le_star/2+La_star+(zc_star-La_star).*(La_star+2*Lc_star-zc_star)/...
    (2*Lc_star)) - rho_v*g*sin(phi(ip))*(Le_star+zc_star);

p_le_star_A = p_sat - delta_p_vt_A + ...
    mu_l*qc_CaseA_g(ip)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (ze_star.^2/(2*Le_star) + ze_star - Lc_star/2 - La_star) + ...
    rho_l*g*sin(phi(ip))*(La_star + Lc_star - ze_star);
p_la_star_A = p_sat - delta_p_vt_A + ...
    mu_l*qc_CaseA_g(ip)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (za_star - Lc_star/2 - La_star) + ...
    rho_l*g*sin(phi(ip))*(La_star + Lc_star - za_star);
p_lc_star_A = p_sat - delta_p_vt_A - ...
    mu_l*qc_CaseA_g(ip)/(K(ig)*rho_l*Aw_star(ig)*h_lv)*...
    (La_star + Lc_star - zc_star).^2/(2*Lc_star) + ...
    rho_l*g*sin(phi(ip))*(La_star + Lc_star - zc_star);

z_star = [ze_star za_star zc_star];
p_v_star_CaseA = [p_ve_star_A p_va_star_A p_vc_star_A];
p_l_star_CaseA = [p_le_star_A p_la_star_A p_lc_star_A];
end 

for i=i  % Analytical  Case B
[K1_star,K2_star,K4_star,Ca_star] = A_FullMassFlowRate_CaseAB(mu,tp,tg,RL_star,RL,Rg_star(ig),RgArr(ig),P_star,tXi0A(ig,1),tXi0A(ig,2),N_An);

delta_p_vt_B = 8*mu_v*qc_CaseB_g(ip)*Leff_star/(pi*rho_v*RL_star^4*h_lv)+...
    rho_v*g*sin(phi(ip))*Lt_star; %is this right???

p_ve_star_B = p_sat - 8*mu_v*qc_CaseB_g(ip)/(pi*rho_v*RL_star^4*h_lv)*...
    (ze_star.^2/(2*Le_star)+ze_star+Le_star/2) - rho_v*g*sin(phi(ip))*...
    (Le_star+ze_star);
p_va_star_B = p_sat - 8*mu_v*qc_CaseB_g(ip)/(pi*rho_v*RL_star^4*h_lv)*...
    (za_star+Le_star/2) - rho_v*g*sin(phi(ip))*(Le_star+za_star);
p_vc_star_B = p_sat - 8*mu_v*qc_CaseB_g(ip)/(pi*rho_v*RL_star^4*h_lv)*...
    (Le_star/2+La_star+(zc_star-La_star).*(La_star+2*Lc_star-zc_star)/...
    (2*Lc_star)) - rho_v*g*sin(phi(ip))*(Le_star+zc_star);

p_le_star_B = p_sat - delta_p_vt_B + ...
    (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv))...
    *mu_l*qc_CaseB_g(ip)/(K1_star + K2_star - Ca_star*K4_star).*...
    (ze_star.^2./(2*Le_star) + ze_star - Lc_star/2 - La_star)+...
    rho_l*g*sin(phi(ip))*(La_star + Lc_star - ze_star);
p_la_star_B = p_sat - delta_p_vt_B + ...
    (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv))...
    *mu_l*qc_CaseB_g(ip)/(K1_star + K2_star - Ca_star*K4_star)*...
    (za_star - Lc_star/2 - La_star)+...
    rho_l*g*sin(phi(ip))*(La_star + Lc_star - za_star);
p_lc_star_B = p_sat - delta_p_vt_B - ...
    (4*mu*K4_star/(pi*rho_v*RL_star^3*h_lv) - 1/(2*Ng*rho_l*h_lv))...
    *mu_l*qc_CaseB_g(ip)/(K1_star + K2_star - Ca_star*K4_star)*...
    (La_star + Lc_star - zc_star).^2/(2*Lc_star)+...
    rho_l*g*sin(phi(ip))*(La_star + Lc_star - zc_star);

p_v_star_CaseB = [p_ve_star_B p_va_star_B p_vc_star_B];
p_l_star_CaseB = [p_le_star_B p_la_star_B p_lc_star_B];

end

% timestr = datetime('2025-03-22','Format','yyyy-MM-dd');
% filename = sprintf('pzplot_T%i_Rgs%1.2i_La_star%1.2i_%s.mat',T,Rg_star(ig)*1E3,timestr)
% save(['Pressure_Plots\' filename]) % Saving Xi Numerical Results

hfig = figure('DefaultAxesFontSize',12); % save the figure handle in a variable
fsLgd = 14
fsAxi = 14*1.5
lw = 1.4 

hold on; grid on; box on
plot(z_star  ,p_v_starA*0.001,'-',Color='#D95319'    ,LineWidth = lw,DisplayName='C$_0$');
plot(z_star  ,p_l_starA*0.001,'-',Color='#D95319'    ,LineWidth = lw,HandleVisibility='off');
plot(z_star,p_v_star_CaseA*0.001,'--',Color='#77AC30',LineWidth = lw,DisplayName='A');
plot(z_star,p_l_star_CaseA*0.001,'--',Color='#77AC30',LineWidth = lw,HandleVisibility='off');
plot(z_star,p_v_star_CaseB*0.001,'-.',Color='black',LineWidth = lw,DisplayName='B');
plot(z_star,p_l_star_CaseB*0.001,'-.',Color='black',LineWidth = lw,HandleVisibility='off');
% xline(0,'-',{'Adiabatic','begin'},LabelVerticalAlignment='bottom');
% xline(La_star(iz),HandleVisibility='off')
xlabel('$z^*$ (m)'); ylabel('$p^*$ (kPa)'); 
xlim([-0.1 0.2])

    titleCell{1,1} ='$T_{\rm{o}}^*=60^{\circ}$C, $L_{\rm{a}}^*=0.1$~m, $\phi = \pi/2$';
    lgd = legend(Location='best',Orientation='horizontal',Box='off',FontSize=fsLgd)    
    title(lgd,titleCell)
    

hfigname = ['Effects_Gravity\GravityEff_p_star_vs_z_star_v2' sprintf('_T%i',T) '_Rgstr' num2str(Rg_star(ig)*100000)];

picturewidth = 17; % set this parameter and keep it forever
hw_ratio = 0.5; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hfig,hfigname,'-dpdf','-painters')
% print(hfig,hfigname,'-dpng','-painters')





%% FUNCTION SECTION
function [tXi00,tXi01] = XiRootAndrew(rho,mu,tp,tg,RL,Rg,N,tXi00_ig)

tXi00 = fzero(@(tXi00) LeadingOrderFun(rho,mu,tp,tg,RL,Rg,N,tXi00),tXi00_ig);
tXi01 = fzero(@(tXi01)   FirstOrderFun(rho,mu,tp,tg,RL,Rg,N,tXi00,tXi01), 1);

    function f = LeadingOrderFun(rho,mu,tp,tg,RL,Rg,N,tXi00)
        [Qv00,Ql00,Qv01,Ql01] = A_FullMassFlowRate(mu,tp,tg,RL,Rg,tXi00,1,N);
        f = tXi00 + rho/mu*Qv00/Ql00;
    end

    function f = FirstOrderFun(rho,mu,tp,tg,RL,Rg,N,tXi00,tXi01)
        [Qv00,Ql00,Qv01,Ql01] = A_FullMassFlowRate(mu,tp,tg,RL,Rg,tXi00,tXi01,N);
        f = tXi01 - rho/(mu*Ql00)*(Ql01*Qv00/Ql00 - Qv01);
    end

end

function Model = ModelingFunction(tp,tg,RL,Rg,Rm0,xm)
% Creat the model and mesh for the XiRoot / Qv0 Solver Function 
Model = createpde(1);  % model start - create the pde model

dl = zeros(10,7); % initialize & creat decomposed geometry matrix (DGM)
% Row #    1  2             3             4             5             6  7  8   9  10
dl(:,1) = [2; RL*cos(tp/2); 0;            RL*sin(tp/2); 0;            1; 0; 0;  0; 0  ];
dl(:,2) = [2; 0;            xm+Rm0;       0;            0;            1; 0; 0;  0; 0  ];
dl(:,3) = [2; xm+Rm0;       Rg;           0;            0;            2; 0; 0;  0; 0  ];
dl(:,4) = [2; RL*cos(tg/2); Rg*cos(tg/2); RL*sin(tg/2); Rg*sin(tg/2); 0; 2; 0;  0; 0  ];
dl(:,5) = [1; Rg;           Rg*cos(tg/2); 0;            Rg*sin(tg/2); 2; 0; 0;  0; Rg ];
dl(:,6) = [1; RL*cos(tg/2); RL*cos(tp/2); RL*sin(tg/2); RL*sin(tp/2); 1; 0; 0;  0; RL ];
dl(:,7) = [1; xm+Rm0;       RL*cos(tg/2); 0;            RL*sin(tg/2); 1; 2; xm; 0; Rm0];
geometryFromEdges(Model,dl); % creats model geometry

end

function SolverOutput = XiNumericalSolver(Model,rho,mu,tXiArg,Type)
%XiRoot / Qv0 Solver Function
k1 = 1;
s1 = 1;
k2 = 1/mu;
s2 = tXiArg;

specifyCoefficients(Model,"m",0,"d",0,"c",k1,"a",0,"f",s1,"Face",1);
specifyCoefficients(Model,"m",0,"d",0,"c",k2,"a",0,"f",s2,"Face",2);
applyBoundaryCondition(Model,"neumann"  ,"Edge",[1,2,3],"q",0,"g",0);
applyBoundaryCondition(Model,"dirichlet","Edge",[4,5,6],"h",1,"r",0);
results = solvepde(Model);

[A,AE] = area(Model.Mesh); % Obtain area of each element (or get the dS):

Trnglp1 = Model.Mesh.Nodes(:,Model.Mesh.Elements(1,:)); % Obtain elements column-wise index:
Trnglp2 = Model.Mesh.Nodes(:,Model.Mesh.Elements(2,:));
Trnglp3 = Model.Mesh.Nodes(:,Model.Mesh.Elements(3,:));
Trnglc = (Trnglp1+Trnglp2+Trnglp3)./ 3;     % Calculate the x y of each element center:

Trnglc_Soln = interpolateSolution(results,Trnglc)'; % Interpolate the soln value at center of each element
QE = AE.*Trnglc_Soln;
E1_Idx = findElements(Model.Mesh,"region","Face",1);
E2_Idx = findElements(Model.Mesh,"region","Face",2);
Qv_tilde_0 = sum(QE(E1_Idx));
Ql_tilde_0 = sum(QE(E2_Idx));

% Ploting Domain & Mesh
% figure; tiledlayout(2,1);
% nexttile % plot the domain
% pdegplot(Model,"VertexLabels","on","EdgeLabels","on","FaceLabels","on");
% title('Domain'); box on; axis equal; % plot the solution
% xlim([0,Rg]);ylim([0,RL*sin(tp/2)]);
% nexttile % plot the mesh
% pdeplot(Model);
% title('Mesh'); box on; axis equal; % plot the solution
% xlim([0,Rg]);ylim([0,RL*sin(tp/2)]);
% pdeplot(results.Mesh,XYData=results.NodalSolution,Mesh="on")

switch Type
    case 'rootFinding'
        SolverOutput = Qv_tilde_0 + (1/rho) * Ql_tilde_0;
        fprintf('-- err = %d\n',SolverOutput);
    case 'tQvCalculate'
        SolverOutput = Qv_tilde_0;
        fprintf('-- tilde_Qv = %d\n\n',SolverOutput);
    case 'velocityField'
        SolverOutput = results;
end

end

function Sigma_Qv = QvRootNumerical(SC,SL,tm0,Rm0,tXi0N,tQv0N,Qv0Arg,N_pt,Case)

rho = SC.rho;
mu = SC.mu;
Ng = SC.Ng;
RL = SC.RL;
Aw = SC.Aw;
Ca_hat = SC.Ca_hat;

Le = SL.Le;
Lc = SL.Lc;
delta_Bo = SL.d_Bo;

p_ve = 1 - 8*Ng*Qv0Arg*Le/(pi*RL^4) - 1*Ng*rho*Qv0Arg*Le/(mu*Aw)+delta_Bo*Le;
delta_pc = 8*Ng*Qv0Arg*Lc/(pi*RL^4) + 1*Ng*rho*Qv0Arg*Lc/(mu*Aw)-delta_Bo*Lc;

% tsrN = [p_ve delta_pc]

%Case La Only 2024/08/28

switch Case
    case 'C1'
        % % interpolate denser array tm
        Rm0_e = 1/(Ca_hat*p_ve);
        Rm0_e = min(max(Rm0),Rm0_e);
        Rm0_c = 1/(Ca_hat*delta_pc);
        % tsr = [Rm0_e Rm0_c] 
        Rm0_c = min(max(Rm0),Rm0_c);
        tm0_e = interp1(Rm0,tm0,Rm0_e);
        tm0_c = interp1(Rm0,tm0,Rm0_c);        
        tm0q = linspace(tm0_e,tm0_c,N_pt);
        tXi0Nq = interp1(tm0,tXi0N,tm0q);
        Rm0q = interp1(tm0,Rm0,tm0q);
        tQv0Nq = interp1(tm0,tQv0N,tm0q);
        % Qlt0q = interp1(tm0_IA,Qlt0_IA,tm0q);
        delta_p0 = 1./(Rm0q*Ca_hat);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;

    case 'C0'
        % % interpolate denser array tm
        [~,I0] = min(abs(tm0-deg2rad(0)));
        tm0q = linspace(tm0(I0),tm0(I0),N_pt);
        tXi0Nq = linspace(tXi0N(I0),tXi0N(I0),N_pt);
        Rm0q = linspace(Rm0(I0),Rm0(I0),N_pt);
        tQv0Nq = linspace(tQv0N(I0),tQv0N(I0),N_pt);
        delta_p0 = linspace(p_ve,delta_pc,N_pt);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;

    case 'C0_57'
        % % interpolate denser array tm
        [~,I57] = min(abs(tm0-deg2rad(57)));
        tm0q = linspace(tm0(I57),tm0(I57),N_pt);
        tXi0Nq = linspace(tXi0N(I57),tXi0N(I57),N_pt);
        Rm0q = linspace(Rm0(I57),Rm0(I57),N_pt);
        tQv0Nq = linspace(tQv0N(I57),tQv0N(I57),N_pt);
        delta_p0 = linspace(p_ve,delta_pc,N_pt);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;

    case 'C0_intMax'
        p_ve = 1 ;
        delta_pc = 0;
        [~,I0] = min(abs(tm0-deg2rad(0)));
        tm0q = linspace(tm0(I0),tm0(I0),N_pt);
        tXi0Nq = linspace(tXi0N(I0),tXi0N(I0),N_pt);
        Rm0q = linspace(Rm0(I0),Rm0(I0),N_pt);
        tQv0Nq = linspace(tQv0N(I0),tQv0N(I0),N_pt);
        delta_p0 = linspace(p_ve,delta_pc,N_pt);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;
end


end

function Sigma_Qv = QvRootAndrew(SC,SL,tXi00,tXi01,Qv00,Qv01,Qv0Arg,Case)

rho = SC.rho;
mu = SC.mu;
Ng = SC.Ng;
R_l = SC.RL;
t_g = SC.tg;
Aw = SC.Aw;
Ca_hat = SC.Ca_hat;
Le = SL.Le;
Lc = SL.Lc;
delta_Bo = SL.d_Bo;

lower_limit = 1 - 8*Ng*Qv0Arg*Le/(pi*R_l^4) - Ng*rho*Qv0Arg*Le/(mu*Aw) + delta_Bo*Le; % p_ve;
upper_limit = 8*Ng*Qv0Arg*Lc/(pi*R_l^4) + Ng*rho*Qv0Arg*Lc/(mu*Aw) - delta_Bo*Lc; % delta_pc = p_vc - p_lc

% tsrA = [lower_limit upper_limit]

switch Case
    case 'C1'
        t_m0 = @(delta_p0) asin(delta_p0*Ca_hat*R_l*sin(t_g/2))-t_g/2;

    case 'C0'
        t_m0 = @(delta_p0) 0;

    case 'C0_57'
        t_m0 = @(delta_p0) deg2rad(57);
        
    case 'C0_intMax'
        lower_limit = 1; % p_ve;
        upper_limit = 0; % delta_pc;
        t_m0 = @(delta_p0) 0;
end

integrand = @(delta_p0) ones(size(delta_p0))/...
    (delta_Bo+Qv0Arg*(1-tXi00-tXi01*t_m0(delta_p0))/...
    (Qv00+Qv01*t_m0(delta_p0)));

Sigma_Qv = integral(integrand,lower_limit,upper_limit) - 1;

end

function Sigma_Qv = QvRootNumerical_g(SC,SL,tm0,Rm0,tXi0N,tQv0N,Qv0Arg,N_pt,Case)

rho = SC.rho;
mu = SC.mu;
Ng = SC.Ng;
RL = SC.RL;
Aw = SC.Aw;
Ca_hat = SC.Ca_hat;

Le = SL.Le;
Lc = SL.Lc;
delta_Bo = SL.d_Bo;

p_ve = 1 - 8*Ng*Qv0Arg*Le/(pi*RL^4) - 1*Ng*rho*Qv0Arg*Le/(mu*Aw)+delta_Bo*Le;
delta_pc = 8*Ng*Qv0Arg*Lc/(pi*RL^4) + 1*Ng*rho*Qv0Arg*Lc/(mu*Aw)-delta_Bo*Lc;

% tsrN = [p_ve delta_pc]

%Case La Only 2024/08/28

switch Case
    case 'C1'
        % % interpolate denser array tm
        Rm0_e = 1/(Ca_hat*p_ve);
        Rm0_e = min(max(Rm0),Rm0_e);
        if Rm0_e < min(Rm0)
            Rm0_e = min(Rm0);
        end
        Rm0_c = 1/(Ca_hat*delta_pc);
        if Rm0_c > max(Rm0)
            Rm0_c = max(Rm0);
        end
        if Rm0_c < min(Rm0)
            Rm0_c = min(Rm0);
        end
        % tsr = [Rm0_e Rm0_c]
        tm0_e = interp1(Rm0,tm0,Rm0_e);
        tm0_c = interp1(Rm0,tm0,Rm0_c);
        tm0q = linspace(tm0_e,tm0_c,N_pt);
        tXi0Nq = interp1(tm0,tXi0N,tm0q);
        Rm0q = interp1(tm0,Rm0,tm0q);
        tQv0Nq = interp1(tm0,tQv0N,tm0q);
        % Qlt0q = interp1(tm0_IA,Qlt0_IA,tm0q);
        delta_p0 = 1./(Rm0q*Ca_hat);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;

        if isnan(Sigma_Qv)
            1
        end 


    case 'C0'
        % % interpolate denser array tm
        [~,I0] = min(abs(tm0-deg2rad(0)));
        tm0q = linspace(tm0(I0),tm0(I0),N_pt);
        tXi0Nq = linspace(tXi0N(I0),tXi0N(I0),N_pt);
        Rm0q = linspace(Rm0(I0),Rm0(I0),N_pt);
        tQv0Nq = linspace(tQv0N(I0),tQv0N(I0),N_pt);
        delta_p0 = linspace(p_ve,delta_pc,N_pt);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;

    case 'C0_57'
        % % interpolate denser array tm
        [~,I57] = min(abs(tm0-deg2rad(57)));
        tm0q = linspace(tm0(I57),tm0(I57),N_pt);
        tXi0Nq = linspace(tXi0N(I57),tXi0N(I57),N_pt);
        Rm0q = linspace(Rm0(I57),Rm0(I57),N_pt);
        tQv0Nq = linspace(tQv0N(I57),tQv0N(I57),N_pt);
        delta_p0 = linspace(p_ve,delta_pc,N_pt);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;

    case 'C0_intMax'
        p_ve = 1 ;
        delta_pc = 0;
        [~,I0] = min(abs(tm0-deg2rad(0)));
        tm0q = linspace(tm0(I0),tm0(I0),N_pt);
        tXi0Nq = linspace(tXi0N(I0),tXi0N(I0),N_pt);
        Rm0q = linspace(Rm0(I0),Rm0(I0),N_pt);
        tQv0Nq = linspace(tQv0N(I0),tQv0N(I0),N_pt);
        delta_p0 = linspace(p_ve,delta_pc,N_pt);
        % Andrew's Integration Formular
        Sigma_Qv = trapz(delta_p0, 1./(delta_Bo+Qv0Arg.*(1-tXi0Nq)./(-tQv0Nq))) - 1;
end


end


