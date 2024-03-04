clear all;
close all;
clc;
load('ECP_values.mat');
global M sigma
%Physical system parameters
J_1 = ECP_values(1);            % Disk 1 inertia kgm^2
J_2 = ECP_values(2);            % Disk 2 inertia kgm^2
J_3 = ECP_values(3);            % Disk 3 inertia kgm^2
k_1 = ECP_values(4);            % Shaft 1-2 stiffness Nm/rad
k_2 = ECP_values(5);            % Shaft 2-3 stiffness Nm/rad
b_1 = mean(ECP_values([6 7]));  % Disk 1 damping and friction Nms/rad
b_2 = mean(ECP_values([8 9]));  % Disk 2 damping and friction Nms/rad
b_3 = mean(ECP_values([10 11]));% Disk 3 damping and friction Nms/rad
T_Cp = ECP_values(12);          % Disk 1 Coulomb friction in positive direction
T_Cm = ECP_values(13);          % Disk 1 Coulomb friction in negative direction
atan_scale = 100;               % Sign approximation factor
w_th = 0.75;                    % Threshold angular velocity rad/s


%The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;                    % Sampling period
sigma_meas = 0.0093*eye(3);     % Measurements covariance matrix

%% Question 1
% Inputs
syms u1(t) u2(t) % Symbolic input declearation
ST_input(1,:) = {u1,u2}; % symbolic variable
ST_input(2,:) = {'u1\\left(t\\right)','u2\\left(t\\right)'}; % LaTeX expression


% Measurements
syms y1(t) y2(t) y3(t) % Symbolic measurement declearation
ST_meas(1,:) = {y1,y2,y3}; % symbolic variable
ST_meas(2,:) = {'y1\\left(t\\right)','y2\\left(t\\right)','y3\\left(t\\right)'}; % LaTeX expression


% Unknowns
syms theta1(t) omega1(t) theta2(t) omega2(t) theta3(t) omega3(t) dtheta1(t) domega1(t) dtheta2(t) domega2(t) dtheta3(t) domega3(t) d(t) % Symbolic unknowns declearation
ST_unknowns(1,:) = {theta1,omega1,theta2,omega2,theta3,omega3,dtheta1,domega1,dtheta2,domega2,dtheta3,domega3,d}; % symbolic variable
ST_unknowns(2,:) = {...
	'theta1\\left(t\\right)','omega1\\left(t\\right)',...
    'theta2\\left(t\\right)','omega2\\left(t\\right)',...
    'theta3\\left(t\\right)','omega3\\left(t\\right)',...
	'\\dot{theta1}\\left(t\\right)','\\dot{omega1}\\left(t\\right)',...
    '\\dot{theta2}\\left(t\\right)','\\dot{omega2}\\left(t\\right)',...
    '\\dot{theta3}\\left(t\\right)','\\dot{omega3}\\left(t\\right)','d\\left(t\\right)'}; % LaTeX expression


% Parameters
syms J1 J2 J3 k1 k2 b1 b2 b3 % Symbolic parameters declaration
ST_parameters(1,:) = {J1,J2,J3,k1,k2,b1,b2,b3}; % symbolic variable
ST_parameters(2,:) = {'J1','J2','J3','k1','k2','b1','b2','b3'}; % LaTeX expression


% Enter the constraints of the system
% ST_cons(1.:) comprise the Matlab names of constraints
% ST_cons(2.:) comprise the Latex 
% ST_cons(3,:) comprise a cell array list with the expressions of the constraints

ST_cons(1,:) = {'c1','c2','c3','c4','c5','c6','d7','d8','d9','d10','d11','d12','m13','m14','m15'}; % Constraint names
ST_cons(2,:) = {'c_1','c_2','c_3','c_4','c_5','c_6','d_7','d_8','d_9','d_10','d_11','d_12','m_13','m_14','m_15'}; % Constraint latex names
ST_cons(3,:) = {...
    0 == dtheta1 - omega1, ...
    0 == J1*domega1 - u1 + b1*omega1 + k1*(theta1 - theta2) + d, ...
    0 == dtheta2 - omega2, ...
    0 == J2*domega2 - u2 + b2*omega2 + k1*(theta2 - theta1) + k2*(theta2 - theta3), ...
    0 == dtheta3 - omega3,...
    0 == J3*domega3 + b3*omega3 + k2*(theta3 - theta2),...
    0 == dtheta1 - diff(theta1,t),...
    0 == dtheta2 - diff(theta2,t),...
    0 == dtheta3 - diff(theta3,t),...
    0 == domega1 - diff(omega1,t),...
    0 == domega2 - diff(omega2,t),...
    0 == domega3 - diff(omega3,t),...
    0 == y1 - theta1,...
    0 == y2 - theta2,...
    0 == y3 - theta3};
% NOTE that "diff(.,t)" is used as the differential operator, D == d/dt.


ST_canfail=[1:6,13:15]; % constraint d7,d8,d9,d10,d11,d12 cannot fail
ST_domains=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
% NOTE NOTE NOTE:
% Incidence matrix MUST have colums ordered as: 
% input, measurements, unknowns
% an x in the book is a "-1" in the software tool

  
% Automatic incMatrix generation
cons_oneway = {[],[],[],[],[],[],{theta1},{theta2},{theta3},{omega1},{omega2},{omega3},[],[],[]};
ST_IncMat = sa_incidencematrix(ST_cons,...
                                    ST_input,ST_meas,...
                                    ST_unknowns,cons_oneway);
						
								
ST_sys =...
    sa_create(ST_IncMat,ST_cons,...
    ST_input, ST_meas,ST_unknowns,...
    ST_domains, ST_canfail,...
	ST_parameters);



sa_disp(ST_sys);

ST_sys=sa_match(ST_sys,'rank');

% Autogenerate pdf report - requires LaTeX to be installed
sa_report(ST_sys,'threeDiskOscillator_sys','pdf','analytic_expressions',true);
disp('A report has been generated with the following results:')

disp('Obtained matching:');
sa_disp(ST_sys, 't');

disp('Parity relation symbolic form:')
sa_disp(ST_sys, 's');

disp('Parity relation analytic form:')
sa_disp(ST_sys, 'a');
%% State space representation
A = [  0          1             0              0        0        0;
    -k_1/J_1  -b_1/J_1       k_1/J_1           0        0        0;
       0          0             0              1        0        0;
    k_1/J_2       0    -(k_1+k_2)/J_2      -b_2/J_2  k_2/J_2     0;
       0          0             0              0        0        1;
       0          0          k_2/J_3           0    -k_2/J_3 -b_3/J_3];
B = [0    0;
    1/J_1 0;
     0    0;
     0   1/J_2;
     0    0;
     0    0];
C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];
D = zeros(3,2);
E_x = [0 1 0 0 0 0]';
E_y = zeros(3,1);
F_x = [0 0 0 0 0;
       0 0 0 1 0;
       0 0 0 0 0;
       0 0 0 0 1;
       0 0 0 0 0;
       0 0 0 0 0];
F_y = [eye(3) zeros(3,2)];
syms s;

%% Question 2

%From question 1 we found the parity relations and the following residuals
f11=tf(0,1);
f12=tf(-k_2,1);
f13=tf([J_3 b_3 k_2],1);
f14=tf(0,1);
f15=tf(0,1);
f21=tf(-k_1,1);
f22=tf([J_2 b_2 k_1+k_2],1);
f23=tf(-k_2,1);
f24=tf(0,1);
f25=tf(-1,1);
F1=[f11 f12 f13 f14 f15];
F2=[f21 f22 f23 f24 f25];
F=[F1;F2];

% Get H_yu transfer function matrix
H_yu = (C*((s*eye(size(A,2)) - A)^-1)*B + D);
% Get H_yd transfer function matrix
H_yd = (C*((s*eye(size(A,2)) - A)^-1)*E_x + E_y);


% Other possible way: working with the left nullspace, but we decided not to use it 
% F = null(H_tilde')';            % Calculate the left nullspace basis
% F = simplify(F);
% F = vpa(F,2);
% % Obtain the "transfer function" object from the symbolic form
% F_tf = tf('s')*zeros(size(F));
% for i = 1:size(F,1)
%     for j = 1:size(F,2)
%         [nn,dd] = numden(F(i,j));
%         nn = sym2poly(nn);
%         dd = sym2poly(dd);
%         F_tf(i,j) = tf(nn,dd);
%     end
% end

% Filter matrix design
w_0 = 10;               % Cut-off frequency
zeta = 0.7;             % Damping (good choices 0.7 - 0.9)
% Use this for second order lowpass filter
LPF_2 = tf(w_0^2,[1 2*zeta*w_0 w_0^2]);
% Q(s) is just a diagonal with the same filter (LPF_1 or LPF_2) everywhere
Q_f = [LPF_2    0;
       0      LPF_2];
% Apply the filter to F(s)
R_f = minreal(Q_f*F);
% 
% Discrete time
sys_d = c2d(R_f,T_s,'tustin');

%% Question 4
H_yf= simplify((C*((s*eye(size(A,2)) - A)^-1)*F_x + F_y));

% Obtain the "transfer function" object from the symbolic form
H_yf_tf = tf('s')*zeros(size(H_yf));
for i = 1:size(H_yf,1)
    for j = 1:size(H_yf,2)
        [nn,dd] = numden(H_yf(i,j));
        nn = sym2poly(nn);
        dd = sym2poly(dd);
        H_yf_tf(i,j) = tf(nn,dd);
    end
end
V_ry=R_f(:,1:3);
H_rf=minreal(V_ry*H_yf_tf);

%Checking weak & strong detectability of the faults

if rank([H_yd H_yf(:,1)])>rank(H_yd)
    [z,gain1]=zero(F(1,:)*[H_yf_tf(:,1);0;0]);
    [z,gain2]=zero(F(2,:)*[H_yf_tf(:,1);0;0]);
    gain=[gain1 gain2];
    if gain1==0 && gain2==0
        display('fault on sensor 1 is weakly detectable')
    else
        display('fault on sensor 1 is strongly detectable')
    end
else
    display('fault on sensor 1 is not detectable')
end
if rank([H_yd H_yf(:,2)])>rank(H_yd)
    [z,gain1]=zero(F(1,:)*[H_yf_tf(:,2);0;0]);
    [z,gain2]=zero(F(2,:)*[H_yf_tf(:,2);0;0]);
    gain=[gain1 gain2];
    if gain1==0 && gain2==0
        display('fault on sensor 2 is weakly detectable')
    else
        display('fault on sensor 2 is strongly detectable')
    end
else
    display('fault on sensor 2 is not detectable')
end
if rank([H_yd H_yf(:,3)])>rank(H_yd)
    [z,gain1]=zero(F(1,:)*[H_yf_tf(:,3);0;0]);
    [z,gain2]=zero(F(2,:)*[H_yf_tf(:,3);0;0]);
    gain=[gain1 gain2];
    if gain1==0 && gain2==0
        display('fault on sensor 3 is weakly detectable')
    else
        display('fault on sensor 3 is strongly detectable')
    end
else
    display('fault on sensor 3 is not detectable')
end
if rank([H_yd H_yf(:,4)])>rank(H_yd)
    [z,gain1]=zero(F(1,:)*[H_yf_tf(:,4);0;0]);
    [z,gain2]=zero(F(2,:)*[H_yf_tf(:,4);0;0]);
    gain=[gain1 gain2];
    if gain1==0 && gain2==0
        display('fault on actuator 1 is weakly detectable')
    else
        display('fault on actuator 1 is strongly detectable')
    end
else
    display('fault on actuator 1 is not detectable')
end
if rank([H_yd H_yf(:,5)])>rank(H_yd)
    [z,gain1]=zero(F(1,:)*[H_yf_tf(:,5);0;0]);
    [z,gain2]=zero(F(2,:)*[H_yf_tf(:,5);0;0]);
    gain=[gain1 gain2];
    if gain1==0 && gain2==0
        display('fault on actuator 2 is weakly detectable')
    else
        display('fault on actuator 2 is strongly detectable')
    end
else
    display('fault on actuator 2 is not detectable')
end

%Coulomb friction known
E_x1=zeros(6,1);
H_yd_2 = (C*((s*eye(size(A,2)) - A)^-1)*E_x1 + E_y);

%same procedure as before, but using H_yd_2 instead of H_yd
%In this case there is no disturbance and all the faults are at least detectable 



%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
h = 0;                  % Put the threshold from GLR here
M = 0;                  % Window size
% Choose residual that is sensetive to fault on y2

syms zz gg mm zzz kk; 
% zz represents the integral variable
% gg reperesents h
% mm reperesents M


% with a magnitude f2 = -0.025
P_F = 0.0001;
P_M = 0.01;
P_D = 1 - P_M;

sigma = ones(3)*0.0093; % declared earlier


pd_zz = (1/(sqrt(2)*gamma(1/2)))*(zz^(-1/2))*(exp(-zz/2));
p_zz = int(pd_zz,zz,2*gg,Inf)
eq_1 = P_F-p_zz==0;
h = double(vpasolve(eq_1,gg))


mu0 = 0; % mean H0
%mu1 = 0.065; % mean H1 will need to calculate from residual and given fault
mu1=0.065;

%sigma = 9.6172e-04;
sigma = var(out.r1)
% Use residual 1



% Calculate g(k)
% use ua or ub, change with most recent version
% use u_ with length 2500
% lambda = MM*((mu1-mu0)^2)/(sigma);
% pd_kk = 0.5*((kk/lambda)^(-0.25)) ...
%     *besseli(-0.5,sqrt(lambda*kk)) ...   
%     *exp(-0.5*(kk+lambda));
% pd_k = int(pd_kk,kk,2*h,Inf);
% eq_2 = pd_k == 0.99;
% 
% disp("Beep poop")
% M = double(vpasolve(eq_2,MM,55))
for M = 1:1000
    lamda = M*((mu1-mu0)^2)/(sigma)
    prob_detect = 1-ncx2cdf(2*h,1,lamda)
    if prob_detect >= P_D
        M
        break
    end
end
%Calculate g(k)
g = zeros(size(out.r1)); 
for k = M:length(out.r1)
    maxsum=0;
    for j=k-M+1:k
        sumin=(sum(abs(out.r1(j:k))-mu0))^2;
        if sumin>maxsum
             maxsum=sumin;
        end
    end
    g(k)=1/(2*sigma*M)*maxsum;
end

t=[0:T_s:10];
figure(1)
subplot(2,1,1)
plot(t,g);
subplot(2,1,2)
plot(t,out.r1);

%% Question 6
%We need to work with the data collected during the experiment, and perform
%again the operations found in subsection GLR, to determine M (in this case
%M=6) 
%load('test1.mat')
%u1=timeseries(u(:,1));
%u1=timeseries(u(:,2));

% In order to understand which sensor is affected by the fault, we tried to
% design 2 secondary residuals in order to make the sensor faults isolable
f31=tf(-k_1,1);
p1=[J_3 b_3 k_2];
p2=[J_2 b_2 k_1+k_2];
p=conv(p1,p2)
p3=[0 0 0 0 -k_2^2];
p=p+p3;
f32=tf(p,p1);
f33=tf(0,1);
f34=tf(0,1);
f35=tf(-1,1);
F_n=[f31 f32 f33 f34 f35;
     f11 f12 f13 f14 f15];
LPF_3 = tf(5^2,[1 2*zeta*w_0 5^2]);
% Q(s) is just a diagonal with the same filter LPF_2
Q_f_n = [LPF_3    0;
       0      LPF_2];
% Apply the filter to F(s)
R_f_n = minreal(Q_f_n*F_n);
sys_d2 = c2d(R_f_n,T_s,'tustin');

% System shows that both residuals signal a fault, this means that the
% fault happens in the second sensor.

%% Question 7
con_sys = ss(A,[B E_x],C,[D zeros(3,1)]);
con_sys_d = c2d(con_sys,T_s,'tustin');
F_d=con_sys_d.A;
G_d=con_sys_d.B(:,1:2);
% State-feedback LQR design
Q_c = diag([2 0 2 0 2.5 0.0024]);
R_c = diag([10 10]);
[K_c,P_c,E_c] = dlqr(F_d,G_d,Q_c,R_c);
C_d=con_sys_d.C;
t_i=[0:T_s:10];
% Scaling of reference
C_3=[0 0 0 0 1 0];
C_ref = pinv(C_3*((eye(6)-F_d+G_d*K_c)^-1)*G_d*K_c);

% Kalman filter with friction estimation - DO NOT MODIFY
F_aug = [F_d G_d(:,1);zeros(1,6) 1];
G_aug = [G_d;0 0];
C_aug = [C zeros(3,1)];
% Kalman gain
L_aug = dlqe(F_aug,eye(7),C_aug,1e-3*eye(7),deg2rad(0.0056)*eye(3));
L_o = L_aug(1:6,:);
L_d = L_aug(7,:);


%% Question 8 Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
B_f=[B(:,1) zeros(6,1)];
va_eig = log(eig(F_d - G_d*K_c))/T_s;     % Continuous time eigenvalues
M_va = place(A,B_f,va_eig);
A_D = A-B_f*M_va;
N_D = pinv(B_f)*B;
B_D = B-B_f*N_D;
C_D = C;

% Then discretise your VA
G_f=[G_d(:,1) zeros(6,1)];
if (rank(G_f)==rank([G_d G_f]))
    disp('Perfect static matching for actuator fault');
else
    disp('Imperfect static matching for actuator fault');
end
va_eig_d = exp(va_eig*T_s);
M_d = place(F_d,G_f,va_eig_d);
F_D = F_d-G_f*M_d;
N_D_d = pinv(G_f)*G_d;
G_D = G_d-G_f*N_D_d;
C_D_d = C_d;
B_change = [1 0;0 0];

%% Question 9: Virtual actuator for sensor fault
%Augmented matrices
E_x_d=con_sys_d.B(:,3);
Fa=[F_d E_x_d; 0 0 0 0 0 0 1];
Ga=[G_d; 0 0];
C_d_a=[C_d zeros(3,1)];
C_f_a=[C_d_a(1,:); 0 0 0 0 0 0 0; C_d_a(3,:)];

if (rank(C_f_a)==rank([C_d_a;C_f_a]))
    disp('Perfect static matching for sensor fault');
else
    disp('Imperfect static matching for sensor fault');
end
% Check observability of the faulty system
if (rank(obsv(Fa,C_f_a))==7)
    disp('Faulty system is observable');
else
    disp('Faulty system is not observable');
end
% Kalman Filter design
St_imp=diag([1 1 0.5 0.5 1 1 0.1]);
Q=eye(7)*T_s;
R=diag([0.0093^2 0.0093^2 0.0093^2])/T_s;
% Discrete time
[L_V,P,Z,eig_obs] = dlqe(Fa,St_imp,C_f_a,Q,R);
F_V = Fa-L_V*C_f_a;
G_V = Ga;
P_V_d = C_d_a*pinv(C_f_a);
C_V_d = C_d_a-P_V_d*C_f_a;

%% Question 8: Simulation for sensor fault (f_u = 0)
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m_time = 8.5;                 % Sensor fault occurence time
% sim('threeDiskOscillatorRig');

%% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
sim('threeDiskOscillatorRig');
sim('threeDiskOscillatorRig_solution');

%% Plot settings
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);

