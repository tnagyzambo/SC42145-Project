%% Setup

clc;
clear;

load("DATA.mat")

G = tf(FWT);

%% Question 3.2

s = tf('s');

G = tf(FWT(:, 1:2));

Wi_1 = ((1 / (16 * pi)) * s + 0.3) / ((1 / (64 * pi)) * s + 1);
Wi_2 = ((1 / (16 * pi)) * s + 0.3) / ((1 / (64 * pi)) * s + 1);
Wi = [Wi_1, 0; ...
      0,    Wi_2];

Wo_1 = (0.05 * s + 0.2) / (0.01 * s + 1);
Wo_2 = (0.05 * s + 0.2) / (0.01 * s + 1);
Wo = [Wo_1, 0;
      0     Wo_2];

Wp = [(s + 5.655) / (3 * s + 0.0005655), 0; ...
      0                                , 0.2];

Wu = [0.01,  0; ...
      0,    (5e-3 * s^2 + 7e-4 * s + 5e-5) / (s^2 + 14e-4 * s + 10^-6)];

P = [zeros(2), zeros(2), zeros(2), Wi;     ...
     G,        zeros(2), zeros(2), G;      ...
     Wp * G,   Wp,       Wp,       Wp * G; ...
     zeros(2), zeros(2), zeros(2), Wu;     ...
    -G,       -eye(2),  -eye(2),  -G];

P.InputName = {'Beta_delta_i (deg)', 'tau_e_delta_i (Nm)', 'Beta_delta_o (deg)', 'tau_e_delta_o (Nm)', 'w_1', 'w_2', 'Beta (deg)', 'tau_e (Nm)'};
P.OutputName = {'Omega_delta_i (rad/s)', 'z_delta_i (m)', 'Omega_delta_o (rad/s)', 'z_delta_o (m)', 'z1_1', 'z1_2', 'z2_1', 'z2_2', 'v(1)', 'v(2)'};

%% Question 3.3

figure(1)
bodemag(Wi_1)
%exportgraphics(gcf, 'images/SC42145_q3_3_wi.png', 'Resolution', 600)

figure(2)
bodemag(Wo_1)
%exportgraphics(gcf, 'images/SC42145_q3_3_wo.png', 'Resolution', 600)

H = ultidyn("H", 1);

P_reduced = P([9, 10], [7, 8]);

figure(3)
bodemag(P_reduced * H); % THIS IS WRONG APPLY H AS PER BLOCK DIAGRAM
%exportgraphics(gcf, 'images/SC42145_q3_3_response.png', 'Resolution', 600)

figure(4)
sigma(P_reduced * H); % THIS IS WRONG APPLY H AS PER BLOCK DIAGRAM
%exportgraphics(gcf, 'images/SC42145_q3_3_sigma.png', 'Resolution', 600)

%% Question 3.5

load('MIXED_SENS_CONT.mat');

P.InputName = {'delta_u(1)', 'delta_u(2)', 'delta_u(3)', 'delta_u(4)', 'w(1)', 'w(2)', 'u(1)', 'u(2)'};
P.OutputName = {'delta_y(1)', 'delta_y(2)', 'delta_y(3)', 'delta_y(4)', 'z1(1)', 'z1(2)', 'z2(1)', 'z2(2)', 'v(1)', 'v(2)'};

Delta = ultidyn("H", 4);
Delta.InputName = {'delta_y(1)', 'delta_y(2)', 'delta_y(3)', 'delta_y(4)'};
Delta.OutputName = {'delta_u(1)', 'delta_u(2)', 'delta_u(3)', 'delta_u(4)'};

P_1 = connect(P, K_2, Delta, {'w(1)', 'w(2)'}, {'z1(1)', 'z1(2)', 'z2(1)', 'z2(2)'});
P_1 = minreal(P_1);

N = lft(P, K_2);
M = lft(Delta, N);