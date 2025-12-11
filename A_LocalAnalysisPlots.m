function [theta,Vzv,dVzvdr] = A_LocalAnalysisPlots(mu,t_p,t_g,R_l,R_g,P_star,h_lv,rho_v,rho_l,xi_tilde,xi_tilde1,N,Ng,tm0,Qv00,Ql00,Qv01,qt_star_c0,qt_star_c1,Case)

R = R_l/R_g; %ratio of radii

[l,n,m,alpha,gamma,beta,Bm_hat,Al_hat,Cn_hat,Dl_hat,En_hat] = A_LocalAnalysisCalc(mu,t_p,t_g,R_l,R_g,xi_tilde,xi_tilde1,N);


%change theta_v and theta_l for global problem
% theta = linspace(0,t_p/2,5000);
theta2 = t_g/2-logspace(-2,-8,2000);
theta3 = t_g/2+logspace(-8,-4,1000);
theta1 = linspace(0,min(theta2),1000) ;
theta4 = linspace(max(theta3),t_p/2,1000) ;
theta = [theta1 theta2 theta3 theta4];

eta = (cos(theta)-cos(t_g/2))/sin(t_g/2); % Eqn 4.20


factorV = qt_star_c0/(2*Ng*h_lv*rho_v*Qv00*P_star^2)

factorL = qt_star_c0/(2*Ng*h_lv*rho_l*Ql00*P_star^2)
factorL2V = factorL * (mu*xi_tilde1)

v_zl_starC0 = -qt_star_c0/(2*Ng*h_lv*rho_l*Ql00*P_star^2)*... 
    (sum(Bm_hat.*(-1).^m.*((1+exp(-2*theta'*beta))./...
    (exp(beta.*(t_g/2-theta'))+exp(-beta.*(t_g/2+theta')))),2)' +...
    sum(Cn_hat.*(R.^(2*gamma)-1).*cos(theta'*gamma),2)');

v_zv_starC0 = qt_star_c0/(2*Ng*h_lv*rho_v*Qv00*P_star^2)*... % Eqn 6.68
    sum(Al_hat.*cos(theta'*alpha),2)';


V_zv00 = sum(Al_hat.*cos(theta'*alpha),2)'; % Eqn 4.49
V_zl00 = V_zv00/(mu*xi_tilde); % Eqn 4.35

% V_zl00_2 = -(sum(Bm_hat.*(-1).^m.*((1+exp(-2*theta_l'*beta))./...
%     (exp(beta.*(t_g/2-theta_l'))+exp(-beta.*(t_g/2+theta_l')))),2)' +...
%     sum(Cn_hat.*(R.^(2*gamma)-1).*cos(theta_l'*gamma),2)');

dV_zv00_dr = R_l/2 + sum(alpha.*Al_hat.*cos(theta'*alpha),2)'/R_l;
dV_zl00_dr = dV_zv00_dr/xi_tilde;

V_zv01 = sum(Dl_hat.*cos(theta'*alpha),2)';  % Eqn 4.74
V_zl01 = sum(En_hat.*(R.^(2*gamma)-1).*cos(theta'*gamma),2)'; % Eqn 4.75

% LHS = R_l*eta.*dV_zv00_dr + V_zv01;
% RHS = mu*(R_l*eta.*dV_zv00_dr + xi_tilde*V_zl01 + xi_tilde1*V_zl00);
% 
% v_zv_starC1 = qt_star_c1/(2*Ng*h_lv*rho_v*(Qv00+tm0.*Qv01)*P_star^2) *...
%     (V_zv00 + tm0.*V_zv01);
% 
% correction =  ( sin(tm0)*cos(theta)+ sqrt(sin(t_g/2)^2-sin(tm0^2)*sin(theta).^2))/(sin(tm0+t_g/2))*R_l - R_l;  
% 
% v_zv_starC1_c = qt_star_c1/(2*Ng*h_lv*rho_v*(Qv00+tm0.*Qv01)*P_star^2) *...
%     (V_zv00 + correction.* dV_zv00_dr + tm0.*V_zv01);

switch Case
    case 'C0'
        Vzv = V_zv00;
        dVzvdr = dV_zv00_dr;

    case 'C1'
        Vzv = V_zv00 + tm0*V_zv01;
        vzv = v_zv_starC1;

    case 'C1c'
        Vzv = V_zv00 + tm0*V_zv01;
        vzv = v_zv_starC1;

    case 'tau'
        Vzv = V_zv00 + tm0*V_zv01;
        vzv = v_zv_starC1;


end 



