%% Initialization
% Clear and close
clc; clear all; close all;

% Load .mat file
load('ECP_values.mat');
load('exp1res.mat');

% Physical system parameters
J1 = ECP_values(1);            % Disk 1 inertia kgm^2
J2 = ECP_values(2);            % Disk 2 inertia kgm^2
J3 = ECP_values(3);            % Disk 3 inertia kgm^2
k1 = ECP_values(4);            % Shaft 1-2 stiffness Nm/rad
k2 = ECP_values(5);            % Shaft 2-3 stiffness Nm/rad
b1 = mean(ECP_values([6 7]));  % Disk 1 damping and friction Nms/rad
b2 = mean(ECP_values([8 9]));  % Disk 2 damping and friction Nms/rad
b3 = mean(ECP_values([10 11]));% Disk 3 damping and friction Nms/rad

% The system states are [theta_1;omega_1;theta_2;omega_2;theta_3;omega_3]
x_0 = [0;0;0;0;0;0];            % Initial conditions
T_s = 0.004;                    % Sampling period
sigma_meas = 0.0093*eye(3);     % Measurements covariance matrix

%% Find the residual generator
syms s U2(s) Y1(s) Y2(s) Y3(s)

% Parameters of the filter
xi = 1.1;       % Damping >\sqrt(2)/2
omega_n = 100;  % Natural frequency of the filter

% 2-order LPF
H(s) = omega_n^2 / (s^2 + 2*xi*omega_n*s + omega_n^2); 

% Parity relations in Laplace domain
eqn1 = U2(s) - Y2(s) * (b2*s + k1 + k2 + s^2*J2) + k1*Y1(s) + k2*Y3(s);
eqn2 = k2*Y2(s) - Y3(s) *(k2 + b3*s + s^2*J3);

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

den = cell2mat(HY3_N1_D.denominator);

% Equation 2
eqn2_f_min = minreal(eqn2_f, [], false);
[N2,D2] = numden(eqn2_f_min);
poly2 = coeffs(N2, [U2(s), Y1(s), Y2(s), Y3(s)]);
poly_den2 = sym2poly(D2);

coefficients_2 = [];
for i=1:length(poly2)
    coefficients_2 = [coefficients_2; poly2(i)];
end

polynomial_Y2_N2 = sym2poly(coefficients_2(1));
polynomial_Y3_N2 = sym2poly(coefficients_2(2));

% Discretization of equation 2
HY2_N2 = tf(polynomial_Y2_N2, poly_den2);
HY3_N2 = tf(polynomial_Y3_N2, poly_den2);

HY2_N2_D = c2d(HY2_N2, T_s, 'tustin');
HY3_N2_D = c2d(HY3_N2, T_s, 'tustin');

num_HY2_N2_D = cell2mat(HY2_N2_D.numerator);
num_HY3_N2_D = cell2mat(HY3_N2_D.numerator);

den2 = cell2mat(HY3_N2_D.denominator);

open_system('Residual_lab.slx')

