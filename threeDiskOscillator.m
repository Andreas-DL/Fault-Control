clear all;
close all;
clc;
load('ECP_values.mat');
% Physical system parameters
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

% The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;                    % Sampling period
sigma_meas = 0.0093*eye(3);     % Measurements covariance matrix

%% STRUCTURAL ANALYSIS
% Define the input
syms u1(t) u2(t);
ST_input(1,:) = {u1, u2};
ST_input(2,:)    = {'u_1', 'u_2'};  % LaTeX % Also as latex(V_d)

% Define the outputs
syms y1(t) y2(t) y3(t);
ST_meas(1,:)    = {y1, y2, y3};
ST_meas(2,:)    = {'y_1', 'y_2', 'y_3'};   % LaTeX

% Define any parameters
syms J1 J2 J3 k1 k2 b1 b2 b3;
ST_par(1,:)    = {J1,J2,J3,k1, k2,b1, b2, b3};
ST_par(2,:)    = {'J1','J2','J3','k1','k2','b1','b2', 'b3'};

% Define the unknown variables
syms theta1(t) dtheta1(t) omega1(t) domega1(t) theta2(t) dtheta2(t) omega2(t) domega2(t) theta3(t) dtheta3(t) omega3(t) domega3(t) d(t);
ST_unknowns(1,:) = {theta1, dtheta1,omega1,domega1, theta2, dtheta2, omega2, domega2, theta3,dtheta3,omega3, domega3, d};
ST_unknowns(2,:) = {'\\theta_1','\\dot{\\theta}_1', '\\omega_1','\\dot{\\omega}_1','\\theta_2','\\dot{\\theta}_2','\\omega_2',...
                    '\\dot{\\omega}_2','\\theta_3','\\dot{\\theta}_3','\\omega_3', '\\dot{\\omega}_3','d'};    % LaTeX

ST_cons(1,:) = {'c1','c2','c3','c4','c5','c6','d7','d8','d9', 'd10',...
    'd11', 'd12', 'm13', 'm14', 'm15'}; % names to be used in the equations
ST_cons(2,:) = {'c_1','c_2','c_3','c_4','c_5','c_6','d_7','d_8','d_9', 'd_{10}',...
    'd_{11}', 'd_{12}', 'm_{13}', 'm_{14}', 'm_{15}'};     % LaTeX notation  
ST_cons(3,:) = {...
    0 == dtheta1 - omega1, ...
    0 == J1*domega1 - u1 + b1*omega1 + k1*(theta1 - theta2) + d, ...
    0 == dtheta2 - omega2, ...
    0 == J2*domega2 - u2 + b2*omega2 + k1*(theta2 - theta1) + k2*(theta2-theta3), ...
    0 == dtheta3 - omega3,...
    0 == J3*domega3 + b3*omega3 + k2*(theta3 - theta2),...
    0 == dtheta1 - diff(theta1,t),...
    0 == domega1 - diff(omega1,t),...
    0 == dtheta2 - diff(theta2,t), ...
    0 == domega2 - diff(omega2,t),...
    0 == dtheta3 - diff(theta3,t), ...
    0 == domega3 - diff(omega3,t),...
    0 == y1 - theta1,...
    0 == y2 - theta2,...
    0 == y3 - theta3};

ST_canfail  = [1 2 3 4 5 6 13 14 15]; % constraints d7 to d12 can't fail
ST_domains  = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Number of constraints

% Constraints that are non-invertible (write the corresponding "known" variable)
ST_cons_oneway = {[],[],[],[],[],[],{theta1},{theta2},{omega1},{omega2},{theta3},{omega3},[],[],[]};
% Compute incidence matrix
ST_IncidenceMatrix = sa_incidencematrix(ST_cons,ST_input,ST_meas,ST_unknowns,ST_cons_oneway);

% Create the system object
threeDiskOscillator_sys = sa_create(ST_IncidenceMatrix,ST_cons,ST_input, ST_meas,ST_unknowns, ST_domains,...
           ST_canfail, ST_par);
% Display the system
sa_disp(threeDiskOscillator_sys);

% Perform a matching
threeDiskOscillator_sys = sa_match(threeDiskOscillator_sys,'rank');

% Autogenerate pdf report - requires LaTeX to be installed
sa_report(threeDiskOscillator_sys,'threeDiskOscillator_sys','pdf',...
            'analytic_expressions',true);
disp('A report has been generated with the following results:')

% Display the results
disp('Obtained matching:');
sa_disp(threeDiskOscillator_sys, 't');

disp('Parity relation symbolic form:')
sa_disp(threeDiskOscillator_sys, 's');
disp('Parity relation analytic form:')
sa_disp(threeDiskOscillator_sys, 'a');


%% Residual filter design
syms s Y1(s) Y2(s) Y3(s) U2(s) U1(s)

eqn1 = U2(s) - Y2(s) * (b_2*s +k_1 + k_2 + s^2*J_2) + k_1*Y1(s) + k_2*Y3(s);
eqn2 = k_2*Y2(s) - Y3(s) *(k_2 + b_3*s + s^2*J_3);

% Parameters of the filter
xi = 0.8;  % Damping >\sqrt(2)/2
omega_n = 20; %natural frequency of the filter

% 2-order LPF
H(s) = omega_n^2 / (s^2 + 2*xi*omega_n*s + omega_n^2); 

% Filter the residuals
eqn1_f = eqn1 * H(s);
eqn2_f = eqn2 * H(s);

eqn1_f_min = minreal(eqn1_f, [], false);
[N1,D1] = numden(eqn1_f_min);
poly1 = coeffs(N1, [U2(s), Y1(s), Y2(s), Y3(s)]);
poly_den1 = sym2poly(D1);

coefficients = [];
for i=1:length(poly1)
    coefficients = [coefficients; poly1(i)];
end

polynomial_U2_N1 = sym2poly(coefficients(4));
polynomial_Y1_N1 = sym2poly(coefficients(3));
polynomial_Y2_N1 = sym2poly(coefficients(2));
polynomial_Y3_N1 = sym2poly(coefficients(1));

% Discretization of equation 1
HU2_N1 = tf(polynomial_U2_N1, poly_den1);
HY1_N1 = tf(polynomial_Y1_N1, poly_den1);
HY2_N1 = tf(polynomial_Y2_N1, poly_den1);
HY3_N1 = tf(polynomial_Y3_N1, poly_den1);

HU2_N1_D = c2d(HU2_N1, T_s, 'tustin');
HY1_N1_D = c2d(HY1_N1, T_s, 'tustin');
HY2_N1_D = c2d(HY2_N1, T_s, 'tustin');
HY3_N1_D = c2d(HY3_N1, T_s, 'tustin');

num_HU2_N1_D = cell2mat(HU2_N1_D.numerator);
num_HY1_N1_D = cell2mat(HY1_N1_D.numerator);
num_HY2_N1_D = cell2mat(HY2_N1_D.numerator);
num_HY3_N1_D = cell2mat(HY3_N1_D.numerator);

den1 = cell2mat(HY3_N1_D.denominator);

% Equation 2
eqn2_f_min = minreal(eqn2_f, [], false);
[N2,D2] = numden(eqn2_f_min);
poly2 = coeffs(N2, [U2(s), Y1(s), Y2(s), Y3(s)]);

coefficients_2 = [];
for i=1:length(poly2)
    coefficients_2 = [coefficients_2; poly2(i)];
end

polynomial_Y2_N2 = sym2poly(coefficients_2(2));
polynomial_Y3_N2 = sym2poly(coefficients_2(1));

poly_den2 = sym2poly(D2);

% Discretization of equation 2
HY2_N2 = tf(polynomial_Y2_N2, poly_den2);
HY3_N2 = tf(polynomial_Y3_N2, poly_den2);

HY2_N2_D = c2d(HY2_N2, T_s, 'tustin');
HY3_N2_D = c2d(HY3_N2, T_s, 'tustin');

num_HY2_N2_D = cell2mat(HY2_N2_D.numerator);
num_HY3_N2_D = cell2mat(HY3_N2_D.numerator);

den2 = cell2mat(HY3_N2_D.denominator);

%% State space representation
A = [0 1 0 0 0 0;
    -k_1/J_1 -b_1/J_1 k_1/J_1 0 0 0;
    0 0 0 1 0 0;
    k_1/J_2 0 -(k_1+k_2)/J_2 -b_2/J_2 k_2/J_2 0;
    0 0 0 0 0 1;
    0 0 k_2/J_3 0 -k_2/J_3 -b_3/J_3];
B = [0 0;
    1/J_1 0;
    0 0;
    0 1/J_2;
    0 0;
    0 0];
C = [  1 0 0 0 0 0
       0 0 1 0 0 0
       0 0 0 0 1 0 ];
D = zeros(3,2);

E_x = [0; 1; 0; 0; 0; 0];
E_y = zeros(3,1);

F_x = [0 0 0 0 0;...
        0 0 0 1 0;...
        0 0 0 0 0;...
        0 0 0 0 1;...
        0 0 0 0 0;...
        0 0 0 0 0];

F_y = [eye(3) zeros(3,2)];
%%
%%%%%%%% SKIP THAT WHEN SIMULATING IN OPEN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Discrete time
% sys_d = [];
% F_d = sys_d.A;
% G_d = sys_d.B;
% 
% % State-feedback LQR design
% Q_c = diag([2 0 2 0 2.5 0.0024]);
% R_c = diag([10 10]);
% K_c = [];
% 
% % Scaling of reference
% C_ref = [];
% 
% % Kalman filter with friction estimation - DO NOT MODIFY
% F_aug = [F_d G_d(:,1);zeros(1,6) 1];
% G_aug = [G_d;0 0];
% C_aug = [C zeros(3,1)];
% % Kalman gain
% L_aug = dlqe(F_aug,eye(7),C_aug,1e-3*eye(7),sigma_meas(1,1).^2*eye(3));
% L_o = L_aug(1:6,:);
% L_d = L_aug(7,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Residual filter design
syms s;
H_yu = simplify(C*inv((s*eye(6) - A))*B + D);
H_yd = simplify(C*inv((s*eye(6) - A))*E_x + E_y);
H_yf = simplify(C*inv((s*eye(6) - A))*F_x + F_y);

H = [H_yu H_yd; eye(2) zeros(2,1)];
xi = 0.8;  % Damping >\sqrt(2)/2
omega_n = 10; %natural frequency of the filter

Q = omega_n^2 / (s^2 + 2*xi*omega_n*s + omega_n^2);

% 2-order LPF
F = simplify(null(H')');
V_ry =  F(:,1:3);




%% Strong and weak detectability
H_rf = simplify(V_ry * H_yf);

f = {'y1', 'y2', 'y3', 'u1', 'u2'};

for i = 1:length(f)
    % Check if the augmented system rank increases
    if rank([H_yd H_yf(:,i)]) > rank(H_yd)
        disp(['Fault on ', f{i}, ' is weakly detectable']);
        strong_detect = F * [H_yf(:,i); zeros(size(B,2),1)];
        % Check if the strong detection condition is met
        if subs(strong_detect, s, 0) ~= 0
             disp(['Fault on ', f{i}, ' is strongly detectable']);
        end
    else
        disp(['Fault on ', f{i}, ' is not weakly detectable']);
    end
end


%% Fault Isolation
% 
% faluts = [fy1 fy2 fy3 fu1 fu2];
H_binary = [0 1 1 0 0;
            1 1 0 0 1];

f_hat = pinv(H_binary) * [eqn1; eqn2];
%% GLR
f_m = [0;-0.025;0];     % Sensor fault vector (added to [y1;y2;y3])
P_F = 0.0001;
P_D = 0.99;

load('C:\Users\giova\OneDrive\Desktop\DTU\Robust and Fault tolerant control\MATLAB\Assignment 1\Fault-Control-main\Residuals_test\exp1res.mat')

% h = chi2inv(1 - P_F, 1)/2; % define the threashold
% 
% % Covariance matrix of the sensor measurements given above
% load('C:\Users\giova\OneDrive\Desktop\DTU\Robust and Fault tolerant control\MATLAB\Assignment 1\Fault-Control-main\Residuals_test\exp1res.mat')
% 
% time = r(:,1);
% mu_0 = 0; % compute the mean before the fault
% mu_1 =  mean(out.r.signals.values(2175:end,2)', 2); % compute the mean after the fault
% sigma = cov(out.r.signals.values(1:2075,1)); % covariance of the output before the fault
% sigma = sqrt(Q);
% 
% syms M X 
% lambda = M*(mu_1 - mu_0)^2/ sigma^2;
% P_D = 0.99;
% func = 1/2 * (X / lambda)^(-1/4)* exp(-(X + lambda) / 2)* besseli(-1/2, sqrt(lambda*X));
% int_func = int(func, X, [2*h, Inf]);
% eqn = P_D - int_func == 0;
% sol = vpasolve(eqn, M, [0, Inf]); % Find M by solving the equation
% 
% M = sol;

syms X h mm
func = (1/(sqrt(2) * gamma(1/2))) * (X^(-1/2)) * (exp(-X/2));
int_func = int(func, X, 2*h, Inf);
eqn = P_F - int_func == 0;
h = double(vpasolve(eqn, h));


mu_0 = 0;                       % mean H0

mu_1 = 0.065;                   % mean H1

%sigma = 9.6172e-04;
sigma = var(out.r.signals.values(:,2));
% Use residual 2

for M = 1:1000
    lambda = M*(mu_1 - mu_0)^2 / sigma;
    prob_detect = 1 - ncx2cdf(2*h, 1, lambda);
    if prob_detect >= P_D
        M;
        break
    end
end
%Calculate g(k)
g = zeros(size(out.r.signals.values(:,2))); 
for k = M:length(out.r.signals.values(:,2))
    maxsum=0;
    for j = k - M + 1:k
        sumin = (sum(abs(out.r.signals.values(j:k)) - mu_0))^2;
        if sumin > maxsum
             maxsum = sumin;
        end
    end
    g(k) = 1/(2 * sigma * M) * maxsum;
end

figure(1)
subplot(2, 1, 1)
plot(out.r.time, g);
legend('$$g(t)$$', 'Interpreter', 'latex')
xlabel('Time[s]')
ylabel('Test statistics')
grid on;
subplot(2, 1, 2)
plot(out.r.time, out.r.signals.values(:,2));
legend('$$r_2(t)$$', 'Interpreter', 'latex')
xlabel('Time[s]')
ylabel('$$r_2(t)$$', 'Interpreter','latex')
grid on

%% Virtual actuator
% Failure in actuator 2
% Do the desing first in continuous time
va_eig_d = [];  % Discrete time eigenvalues
va_eig = log(va_eig_d)/T_s;     % Continuous time eigenvalues
% Then discretise your VA

B_change = [1 0;0 0];

%% Simulation for sensor fault (f_u = 0)
simTime = 45;                   % Simulation duration in seconds
f_u_time = 25;                  % Actuator fault occurence time
detect_time = f_u_time + 3.75;
f_u = [0;0];                    % Actuator fault vector (added to [u1;u2])
u_fault = 0;                    % Disable VA meachanism
f_m_time = 8.5;                 % Sensor fault occurence time
sim('threeDiskOscillatorRig');

%% Simulation for actuator fault (f_m = 0)
f_u = [0;-0.1];                 % Actuator fault vector (added to [u1;u2])
u_fault = 1;                    % Enable VA meachanism
f_m = [0;0;0];                  % Sensor fault vector (added to [y1;y2;y3])
sim('threeDiskOscillatorRig');

%% Plot settings
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);

