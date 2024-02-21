clear all;
close all;
clc;

load('rig_init.mat');

%% Residual filter design
syms s Y1(s) Y2(s) Y3(s) U2(s)

eqn1 = U2(s) - Y2(s) * (b_2*s +k_1 + k_2 + s^2*J_2) + k_1*Y1(s) + k_2*Y3(s);
eqn2 = k_2*Y2(s) - Y3(s) *(k_2 + b_3*s + s^2*J_3);

% Parameters of the filter
xi = 0.8;  % Damping >\sqrt(2)/2
omega_n = 10; %natural frequency of the filter

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

den = cell2mat(HY3_N1_D.denominator);

% Equation 2
eqn2_f_min = minreal(eqn2_f, [], false);
[N2,D2] = numden(eqn2_f_min);
poly2 = coeffs(N2, [U2(s), Y1(s), Y2(s), Y3(s)]);

coefficients_2 = [];
for i=1:length(poly2)
    coefficients_2 = [coefficients_2; poly2(i)];
end

polynomial_Y2_N2 = sym2poly(coefficients_2(1));
polynomial_Y3_N2 = sym2poly(coefficients_2(2));

poly_den2 = sym2poly(D2);

% Discretization of equation 2
HY2_N2 = tf(polynomial_Y2_N2, poly_den2);
HY3_N2 = tf(polynomial_Y3_N2, poly_den2);

HY2_N2_D = c2d(HY2_N2, T_s, 'tustin');
HY3_N2_D = c2d(HY3_N2, T_s, 'tustin');

num_HY2_N2_D = cell2mat(HY2_N2_D.numerator);
num_HY3_N2_D = cell2mat(HY3_N2_D.numerator);

den2 = cell2mat(HY3_N2_D.denominator);

%% GLR


%% Virtual sensor (Optional)


%% Plots

set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);