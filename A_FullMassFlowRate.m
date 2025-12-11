function [Qv00,Ql00,Qv01,Ql01] = A_FullMassFlowRate(mu,t_p,t_g,R_l,R_g,xi_tilde,xi_tilde1,N)
%% Step 1: Define constants
% format long;
%N is related to number of terms in series

R = R_l/R_g; %ratio of radii

%define indices to generate matrices of functions more easily

l = (0:N-1); %alpha
n = (0:N-1); %gamma
m = (1:2*N); %beta (will be summed, so need not have exactly N terms)

r = (0:N-1).'; %alpha eigenfxn
s = (0:N-1).'; %gamma eigenfxn

alpha = 2*pi*l/t_p; %generates alpha values
gamma = pi*(1+2*n)/t_g; %generates gamma values

beta = ((m - 1/2)*pi / log(1/R));

Bm_hat = (4*R_g^2*log(R)*beta-(2*R_l^2*log(R)*beta.^2+...
    (R_g^2-R_l^2)*(4+beta.^2)).*(-1).^m)./...
    (2*beta.^2.*(4+beta.^2)*(log(R))^2);

%% Step 2: Create Function Tables

d = (t_p/(4*pi))*((sin(pi*(l-r)*t_g/t_p)./...
    (l-r))+(sin(pi*(l+r)*t_g/t_p)./(l+r))); %l!=r
d(1:N+1:end) = t_g/4 + (t_p./(8*pi*l)).*sin(2*pi*l*t_g/t_p); %l==r!=0
d(1,1) = t_g/2; %l==r==0

e = (-1).^n.*(1+2*n)*t_g*t_p^2.*cos(r*pi*t_g/t_p)./...
    (pi*(t_p^2*(4*n.^2+4*n+1)-4*r.^2*t_g^2));

f = (2*t_p*log(R)./((1-2*m).^2*pi*t_p^2+16*pi*r.^2*log(R)^2)).*...
    (4*r*log(R).*sin(r*pi*t_g/t_p)-(2*m-1)*t_p.*cos(r*pi*t_g/t_p).*...
    tanh(beta*t_g/2));

g = (-1).^s.*(1+2*s)*t_g*(t_p^2).*cos(l*pi*t_g/t_p)./...
    (pi*(t_p^2*(4*s.^2+4*s+1)-4*l.^2*t_g^2));

h = t_g/4*eye(N);

j = (((t_g/pi)*(-1).^s)./(1+2*s));

k = (r*t_p.*sin(r*pi*t_g/t_p).*cos(l*pi*t_g/t_p)-l*t_p.*...
    cos(r*pi*t_g/t_p).*sin(l*pi*t_g/t_p))./(2*pi*(l.^2-r.^2)); %l!=r
k(1:N+1:end) = (t_p-t_g)/4 - (t_p./(8*pi*l)).*sin(2*pi*l*t_g/t_p); %l==r!=0
k(1,1) = (t_p-t_g)/2; %l==r==0

o = (t_p^2)*(t_p*cot(t_g/2)*sin(pi*t_g/t_p.*r)-2*pi*r.*...
    cos(pi*t_g/t_p.*r))./(8*pi^3*r.^3-2*pi*t_p^2*r);
o(1)=1-cot(t_g/2)*t_g/2;

p = t_p/4*(csc(t_g/2)*(sin(t_g/2+pi*t_g/t_p*(l+r))./(2*pi*(l+r)+t_p)...
    -sin(t_g/2-pi*t_g/t_p*(l+r))./(2*pi*(l+r)-t_p)...
    +sin(t_g/2+pi*t_g/t_p*(l-r))./(2*pi*(l-r)+t_p)...
    -sin(t_g/2-pi*t_g/t_p*(l-r))./(2*pi*(l-r)-t_p))...
    -cot(t_g/2)/pi*(sin(pi*t_g/t_p*(l-r))./(l-r)...
    +sin(pi*t_g/t_p*(l+r))./(l+r)));
p(1:N+1:end) = (64*pi^3*l.^3-4*pi*t_p^2*(1+cos(2*pi*t_g/t_p*l)).*l...
    +cot(t_g/2)*(-32*pi^3*t_g*l.^3+2*pi*t_g*t_p^2*l+t_p^3 ...
    *sin(2*pi*t_g/t_p*l)))./(128*pi^3*l.^3-8*pi*t_p^2*l);
p(1,1) = 1-cot(t_g/2)*t_g/2;

q = t_g*t_p^2/(2*pi)*(1+2*n).*(-1).^n.*(2*cot(t_g/2)*cos(pi*t_g/t_p*r)...
    ./(4*t_g^2*r.^2-(t_p*(1+2*n)).^2)+pi^2*csc(t_g/2)...
    *(-cos(t_g*(1/2-pi/t_p*r))./((2*pi*t_g*r-t_p*(t_g+pi*(1+2*n)))...
    .*(2*pi*t_g*r+t_p*(-t_g+pi*(1+2*n))))...
    -cos(t_g*(1/2+pi/t_p*r))./((2*pi*t_g*r-t_p*(-t_g+pi*(1+2*n)))...
    .*(2*pi*t_g*r+t_p*(t_g+pi*(1+2*n))))));

t = (t_g^3*cot(t_g/2)*(-1).^s)...
    ./(pi*(1+2*s).*(pi*(1+2*s)-t_g).*(pi*(1+2*s)+t_g));

u = t_g*t_p^2/(2*pi)*(1+2*s).*(-1).^s.*(2*cot(t_g/2)*cos(pi*t_g/t_p*l)...
    ./(4*t_g^2*l.^2-(t_p*(1+2*s)).^2)+pi^2*csc(t_g/2)...
    *(-cos(t_g*(1/2-pi/t_p*l))./((2*pi*t_g*l-t_p*(t_g+pi*(1+2*s)))...
    .*(2*pi*t_g*l+t_p*(-t_g+pi*(1+2*s))))...
    -cos(t_g*(1/2+pi/t_p*l))./((2*pi*t_g*l-t_p*(-t_g+pi*(1+2*s)))...
    .*(2*pi*t_g*l+t_p*(t_g+pi*(1+2*s)))))); %exact same as q

v = 16/pi*t_g^3*log(R)^3*(1+2*s).*(-1).^(1+s).*((cot(t_g/2)*log(R)...
    *(3*pi^2*t_g^2*(1-2*m).^2-4*log(R)^2*(pi*(1+2*s)-t_g)...
    .*(pi*(1+2*s)+t_g)))./((t_g^2*(1-2*m).^2+4*log(R)^2*(1+2*s).^2)...
    .*(t_g^2*pi^2*(1-2*m).^2+4*log(R)^2*(pi*(1+2*s)-t_g).^2)...
    .*(t_g^2*pi^2*(1-2*m).^2+4*log(R)^2*(pi*(1+2*s)+t_g).^2))...
    +(pi^3*tanh(beta*t_g/2).*(2*m-1))./(pi^4*t_g^4*(1-2*m).^4 ...
    +8*pi^2*t_g^2*log(R)^2*(1-2*m).^2.*(t_g^2+(pi*(1+2*s)).^2)...
    +16*log(R)^4*(pi*(1+2*s)-t_g).^2.*(pi*(1+2*s)+t_g).^2));

w = -2*pi^2*t_g^2*(-1).^(n+s).*(1+2*n).*(1+2*s)...
    ./(16*pi^4*(n-s).^2.*(1+n+s).^2 ...
    -4*pi^2*t_g^2*(1+2*n.*(1+n)+2*s.*(1+s))+t_g^4);
w(1:N+1:end) = 1/4*(2-2*t_g^2./(t_g^2-4*((pi+2*pi*n).^2))-t_g*cot(t_g/2));

x = log(R)*t_g*(-1).^s.*((pi*t_g*cot(t_g/2)*(1-2*m)...
    -2*log(R)*tanh(beta*t_g/2).*(pi*(1+2*s)-t_g))...
    ./(pi^2*t_g^2*(1-2*m).^2+4*log(R)^2*(pi*(1+2*s)-t_g).^2)...
    -(pi*t_g*cot(t_g/2)*(1-2*m)...
    +2*log(R)*tanh(beta*t_g/2).*(pi*(1+2*s)+t_g))...
    ./(pi^2*t_g^2*(1-2*m).^2+4*log(R)^2*(pi*(1+2*s)+t_g).^2));

y = t_g/4*(-1).^(n+s).*(1./(2*pi*(n-s)+t_g)+1./(2*pi*(n-s)-t_g)...
    -1./(2*pi*(1+n+s)+t_g)-1./(2*pi*(1+n+s)-t_g));

z = pi/2*t_g*t_p^2*csc(t_g/2)*(1+2*s).*(-1).^s...
    .*(cos(t_g*(1/2-pi/t_p*l))./((2*pi*t_g*l-t_p*(t_g+pi*(1+2*s)))...
    .*(2*pi*t_g*l+t_p*(-t_g+pi*(1+2*s))))...
    -cos(t_g*(1/2+pi/t_p*l))./((2*pi*t_g*l-t_p*(-t_g+pi*(1+2*s)))...
    .*(2*pi*t_g*l+t_p*(t_g+pi*(1+2*s)))));
%% Step 3: Sum functions and concatenate submatrices to create matrix M

M1 = d + k; clear d k  %top left 
M2 = (-1*mu*xi_tilde)*((R.^(2*gamma))-1).*e; %top right
M3 = (alpha/R_l).*g; clear g %bottom left
M4 = (-1*xi_tilde)*(gamma/R_l).*((R.^(2*gamma))+1).*h; %bottom right

M = [M1 M2; M3 M4]; %final 2N x 2N matrix
decomp_M = decomposition(M); clear M

clear M1 M2 M3 M4

%% Step 4: Create RHS

%create temp matrix Ts containing all possible terms on RHS of vel eqn in r
Tr = ((mu*xi_tilde)*Bm_hat.*((-1).^m).*f);

%similarly, create Tr for shear stress in s
Ts = ((R_l/2)*(xi_tilde-1)-xi_tilde*(R_l^2-R_g^2)/(4*R_l*log(R)))*j;

%sum rows of temporary matrix and column vector to create f_total (RHS)
f_r = sum(Tr,2); clear Tr %sum over rows
f_s = Ts;
f_total = [f_r; f_s]; %RHS done

clear f_r f_s
%% Step 5: Solve Matrix Equation

a = decomp_M\f_total;
Al_hat = a(1:N).';
Cn_hat = a(N+1:2*N).';

%% Step 6: Create RHS of First Order Correction (FOC)

%VELOCITY

%start with all N-numbered terms (NxN matrix)
%note use of eye so 1D column vectors are only added once
Tr1a = -R_l^2/2*o.*eye(N)-Al_hat.*alpha.*p...
    +mu*xi_tilde*(R_l^2/2*o.*eye(N)-(R_l^2-R_g^2)/(4*log(R))*o.*eye(N)...
    +Cn_hat.*gamma.*((R.^(2*gamma))+1).*q)...
    +mu*xi_tilde1*(Cn_hat.*((R.^(2*gamma))-1).*e); clear e o p q

%then the M-numbered term (NxM matrix, here M=2N)
Tr1b = mu*xi_tilde1*Bm_hat.*((-1).^m).*f; clear f %FIXED

%sum rows of both Tr1a and Tr1b together
f_r1 = sum(Tr1a,2) + sum(Tr1b,2); clear Tr1a Tr1b

%SHEAR STRESS

%first the N-numbered terms, using eye on 1D column vectors as before
Ts1a = (xi_tilde-1)*R_l/2*t.*eye(N)...
    +xi_tilde*(R_l^2-R_g^2)/(4*R_l*log(R))*t.*eye(N)...
    -1/R_l*Al_hat.*alpha.*(z+(alpha-1).*u)...
    +xi_tilde/R_l*Cn_hat.*gamma.*(((gamma-1).*R.^(2*gamma)-gamma-1).*w...
    +(R.^(2*gamma)-1).*y)...
    +xi_tilde1*((R_l/2-(R_l^2-R_g^2)/(4*R_l*log(R)))*j.*eye(N)...
    +1/R_l*Cn_hat.*gamma.*(R.^(2*gamma)-1).*h);
    clear h u w y z
%then the M-numbered terms
Ts1b = -xi_tilde/R_l*(Bm_hat.*beta.*(-1).^m.*(x+beta.*v));
    
f_s1 = sum(Ts1a,2) + sum(Ts1b,2); clear Ts1a Ts1b

%create final forcing

f_total1 = [f_r1; f_s1]; %RHS done

%% Step 7: Solve First-Order Matrix Equation

a1 = decomp_M\f_total1;
Dl_hat = a1(1:N).';
En_hat = a1(N+1:2*N).';


%% Step 8: Resolve Leading Order Qv

Qv00 = (Al_hat(1)*R_l^2*t_p/4)-(R_l^4*t_p/32); %only 1st term @ lead. order

%% Step 9: Resolve Leading Order Ql

K1 = (t_g/32)*(R_l^4-R_g^4-(R_g^2-R_l^2)^2/log(R));

%this had to be rewritten with the new Bm_hat!! The sinh became tanh
K2_terms = -Bm_hat.*... 
    ((2*R_l^2*(-1).^m + beta*R_g^2).*tanh(beta*t_g/2))./...
    (beta.*(4+beta.^2)); %updated to reflect sin and cos simplifications
K2 = sum(K2_terms);

K3_terms = Cn_hat.*sin(t_g/2*gamma).*...
    (R_l^2*R.^(2*gamma).*(gamma-2) - 2*gamma.*R_g^2.*R.^gamma...
    + R_l^2*(2+gamma))./(gamma.*(4-gamma.^2));
K3 = sum(K3_terms);

Ql00 = K1 + K2 + K3;

%% Step 10: Resolve First Order Qv and Ql

Qv1_01 = R_l^2*Dl_hat(1)*t_p/4;

Qv2_01_terms = R_l^2*Al_hat.*(sin(t_g/2*alpha)/tan(t_g/2) - ...
    alpha.*cos(t_g/2*alpha))./(alpha.*(alpha.^2-1));
Qv2_01_terms(1) = R_l^2*Al_hat(1)*(1 - t_g/(2*tan(t_g/2))); %L'Hospital's
Qv2_01 = sum(Qv2_01_terms);

Ql1_01_terms = En_hat.*sin(t_g/2*gamma).*...
    (R_l^2*R.^(2*gamma).*(gamma-2) - 2*gamma.*R_g^2.*R.^gamma...
    + R_l^2*(2+gamma))./(gamma.*(4-gamma.^2));
Ql1_01 = sum(Ql1_01_terms); %same as K3

Ql2_01_terms_m = R_l^2*Bm_hat.*(-1).^m.*(tanh(t_g/2*beta)/tan(t_g/2)-beta)./...
    (beta+beta.^3); %FIXED WRONG DENOMINATOR
Ql2_01_terms_n = R_l^2*Cn_hat.*(R.^(2*gamma)-1).*...
    (gamma.*cos(t_g/2*gamma)-sin(t_g/2*gamma)/tan(t_g/2))./...
    (gamma.*(gamma.^2-1));
Ql2_01 = sum(Ql2_01_terms_m) + sum(Ql2_01_terms_n);

%let Qi01 = Qi1_01 + Qi2_01 since we only need them summed!

Qv01 = Qv1_01 + Qv2_01;
Ql01 = Ql1_01 + Ql2_01;


end


