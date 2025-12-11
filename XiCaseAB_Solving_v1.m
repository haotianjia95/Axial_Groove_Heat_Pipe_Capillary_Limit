clear all; close all; clc; format compact; format  shortE

% Solver Settings ---------------------------------------------------------
N_An = 1000; % number of terms, Andrew's code

Case = 'B'

% Fluid Proporties --------------------------------------------------------
    T = 200;
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


% % theta_m^0 and R_m^0 array ---------------------------------------------

tm0 = [0];

tm0_deno = sin(tm0 + tg/2);

Rm0_star = RL_star*sin(tg/2)./tm0_deno; %  greater than RL*sin(tg/2), no upper limit
xm_star = RL_star*cos(tg/2) - sqrt(Rm0_star.^2 - RL_star^2 * sin(tg/2)^2); % Rm center x

% figure
% plot(tm0,Rm0,'.')
% xline(pi/2 - tg/2,'-','pi/2 - tg/2')
% yline(Rp,':','Rp,min')
% xlabel('\theta_m^0')
% ylabel('R_m^0')

% % Creat Data Structures -------------------------------------------------
timestr = datetime('2025-07-17','Format','yyyy-MM-dd');

switch Case
    case 'A'
        filename = sprintf('XiCaseANumRgs%1.2i-%1.2i_T%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,timestr)
    case 'B'
        filename = sprintf('XiCaseBNumRgs%1.2i-%1.2i_T%i_%s.mat',min(Rg_star)*1E3,max(Rg_star)*1E3,T,timestr)
end


%% Calculate \tilde\xi^00 & \tilde\xi^01 Values
dpldz_Num = NaN(N_geom,length(tm0));

meMax = 0.0025 * P_star ; % max mesh element size

for i_geom = 1:length(Rg_star) % loop counter for different domain shape

    
    dpvdz = 10 % set this as 

    % Analytical Calculation ----------------------------------------------
    % initial guess analytical:
    if i_geom == 1; dpldz_ig = polyval([-0.0032,0.4,-32.85],T) * dpvdz %-48.93 + 0.8744*T - 0.005149*T^2; % -180 for 200C
    else          ; dpldz_ig = 0.8*dpldz_Num(i_geom-1,1);    end % 0.9 for 200C
    % initial guess analytical:

    

    % Numerical Calculation -----------------------------------------------

    fprintf('Geom %i/%i, Iteration %i/%i, Rm0 = %1.2i.\n',i_geom,length(Rg_star),i,length(Rm0_star),Rm0_star);

    % Model = ModelingFunction(tp,tg,RL,Rg,Rm0,xm); % domain only determined by [tp tg RL Rg]


    % Creat the model and mesh for the XiRoot / Qv0 Solver Function
    Model = createpde(1);  % model start - create the pde model

    dl = zeros(10,4); % initialize & creat decomposed geometry matrix (DGM)
    % Row #    1  2                  3                          4                  5                          6  7  8        9  10
    dl(:,1) = [2; xm_star+Rm0_star ; Rg_star(i_geom)          ; 0                ; 0                        ; 1; 0; 0      ; 0; 0       ];
    dl(:,2) = [2; RL_star*cos(tg/2); Rg_star(i_geom)*cos(tg/2); RL_star*sin(tg/2); Rg_star(i_geom)*sin(tg/2); 0; 1; 0      ; 0; 0       ];
    dl(:,3) = [1; Rg_star(i_geom)  ; Rg_star(i_geom)*cos(tg/2); 0                ; Rg_star(i_geom)*sin(tg/2); 1; 0; 0      ; 0; Rg_star ];
    dl(:,4) = [1; xm_star+Rm0_star ; RL_star*cos(tg/2)        ; 0                ; RL_star*sin(tg/2)        ; 0; 1; xm_star; 0; Rm0_star];
    geometryFromEdges(Model,dl); % creats model geometry
    pdegplot(Model,EdgeLabels="on",FaceLabels="on")
    Model.Mesh = generateMesh(Model,'Hmax',meMax); % meshing from DGM

    dpldz_Num(i_geom) = fzero(@(dpldz_arg) XiCaseBNumSolver(Model,RL_star,tp,rho_l,rho_v,mu_l,mu_v,dpvdz,dpldz_arg,Case,'rootFinding'),dpldz_ig);
    % tQv0N(i_geom) = XiCaseBNumSolver(Model,rho,mu,tXi0N(i_geom,i),Case,'tQvCalculate');
    ssQv0N_caseB(i_geom) = XiCaseBNumSolver(Model,RL_star,tp,rho_l,rho_v,mu_l,mu_v,dpvdz,dpldz_Num(i_geom),Case,'tQvCalculate')

end



switch Case
    case 'A'
        tXi00_CaseANum = dpldz_Num/dpvdz
        save(filename,'tXi00_CaseANum') % Saving Xi Numerical Results
    case 'B'
        tXi00_CaseBNum = dpldz_Num/dpvdz
        save(filename,'tXi00_CaseBNum') % Saving Xi Numerical Results
end


%% FUNCTION SECTION
function SolverOutput = XiCaseBNumSolver(Model,RL_star,tp,rho_l,rho_v,mu_l,mu_v,dpvdz,dpldz,CaseArg,Type)
%XiRoot / Qv0 Solver Function

k = mu_l;
s = dpldz;
q = 0.5*RL_star*dpvdz;

specifyCoefficients(Model,"m",0,"d",0,"c",k,"a",0,"f",s,"Face",1);

switch CaseArg
    case 'A'
        applyBoundaryCondition(Model,"dirichlet","Edge",[ 4 ],"h",1,"r",0);
    case 'B'
        applyBoundaryCondition(Model,"neumann"  ,"Edge",[ 4 ],"q",0,"g",q);
end

applyBoundaryCondition(Model,"neumann"  ,"Edge",[ 1 ],"q",0,"g",0);
applyBoundaryCondition(Model,"dirichlet","Edge",[2,3],"h",1,"r",0);
results = solvepde(Model);

[A,AE] = area(Model.Mesh); % Obtain area of each element (or get the dS):

Trnglp1 = Model.Mesh.Nodes(:,Model.Mesh.Elements(1,:)); % Obtain elements column-wise index:
Trnglp2 = Model.Mesh.Nodes(:,Model.Mesh.Elements(2,:));
Trnglp3 = Model.Mesh.Nodes(:,Model.Mesh.Elements(3,:));
Trnglc = (Trnglp1+Trnglp2+Trnglp3)./ 3;     % Calculate the x y of each element center:

Trnglc_Soln = interpolateSolution(results,Trnglc)'; % Interpolate the soln value at center of each element
QE = AE.*Trnglc_Soln;
E1_Idx = findElements(Model.Mesh,"region","Face",1);
Ql = sum(QE(E1_Idx));

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
        Qv = tp/2 * dpvdz * RL_star^4 / (16*mu_v);
        SolverOutput = Ql*rho_l + Qv*rho_v ;
        fprintf('-- err = %d\n',SolverOutput);
    case 'tQvCalculate'
        % Qv = tp/2 * rho_v * dpvdz * RL_star^4 / (16*mu_v);
        SolverOutput = Ql;
        fprintf('-- tilde_Qv = %d\n\n',SolverOutput);
end

end