%% Setup

clc;
clear;

load("DATA.mat")

G = tf(FWT);

%% Question 2.1

syms s

G_11 = poly2sym(cell2mat(G.Numerator(1, 1)), s) / poly2sym(cell2mat(G.Denominator(1, 1)), s);
G_12 = poly2sym(cell2mat(G.Numerator(1, 2)), s) / poly2sym(cell2mat(G.Denominator(1, 2)), s);
G_13 = poly2sym(cell2mat(G.Numerator(1, 3)), s) / poly2sym(cell2mat(G.Denominator(1, 3)), s);
G_21 = poly2sym(cell2mat(G.Numerator(2, 1)), s) / poly2sym(cell2mat(G.Denominator(2, 1)), s);
G_22 = poly2sym(cell2mat(G.Numerator(2, 2)), s) / poly2sym(cell2mat(G.Denominator(2, 2)), s);
G_23 = poly2sym(cell2mat(G.Numerator(2, 3)), s) / poly2sym(cell2mat(G.Denominator(2, 3)), s);
G_sym = [G_11 G_12 G_13; G_21 G_22 G_23];

RGA = G_sym .* pinv(G_sym)';

RGA_1 = double(subs(RGA, s, 0));

RGA_2 = double(subs(RGA, s, 0.3 * 2 * pi));

%% Question 2.2

sys_poles = pole(G);
sys_zeros = tzero(G);

%% Question 2.3

G = tf(FWT(:, 1:2));

s = tf('s');

w_B = 0.3*2*pi;   % Bandwidth
atn = 10e-5; % Attenuation
M = 3;       % H_inf norm bound

Wp11 = ((s / M) + w_B) / (s + (w_B * atn));

Wp = [Wp11, 0; ...
       0,   0.2];

Wu = [0.01,  0,                                                      ; ...
      0,    (5e-3 * s^2 + 7e-4 * s + 5e-5) / (s^2 + 14e-4 * s + 10-6)];
      %0,     0,                                                        0];

Wt = [];

[K, CL, gamma, info] = mixsyn(G, Wp, Wu, Wt);

S = inv(eye(2) + G * K);

w = logspace(-3, 3, 1001) * 2 * pi;
[mag_Si, phase_Si, wout_Si] = bode(S, w);
[mag_Wp, phase_Wp, wout_Wp] = bode(1/Wp, w);
[mag_Wu, phase_Wu, wout_Wu] = bode(Wp*S, w);

% Plotting
f1 = figure(1);
tiled1 = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'none');

% Gain plot
nexttile(tiled1);
p1 = semilogx(w, squeeze(mag_Si(1, 1, :)), 'color', [0.8500 0.3250 0.0980]); % Plant dynamics
hold on
p2 = semilogx(w, squeeze(mag_Wp(1, 1, :)), 'color', [0.9290 0.6940 0.1250]); % Controller dynamics
yline(0, ':', 'color', [0, 0, 0] + 0.25); % Zero gain line
xlim([w(1), w(end)])
ylim([-0.5, 1.5]);
legend('|S|', '|1/Wp|');

nexttile(tiled1);
p1 = semilogx(w, squeeze(mag_Wu(1, 1, :)), 'color', [0.8500 0.3250 0.0980]); % Plant dynamics
hold on
%p2 = semilogx(w, squeeze(mag_Wp(1, 1, :)), 'color', [0.9290 0.6940 0.1250]); % Controller dynamics
yline(0, ':', 'color', [0, 0, 0] + 0.25); % Zero gain line
xlim([w(1), w(end)])
ylim([0, 2]);
legend('|WpS|');
hold off

hinfnorm(S)

figure(2)
sigma(S)

%% Question 2.5

% Plant
G.InputName = {'u(1)'; 'u(2)'};
G.OutputName = 'y';

% Controller
%K.InputName = 'v';
%K.OutputName = 'u';

% Wp
Wp.InputName = 'v';
Wp.OutputName = 'z1';

% Wu
Wu.InputName = 'u';
Wu.OutputName = 'z2';

% Sum block
Sum1 = sumblk('v(1) = w(1) + y(1)');
Sum2 = sumblk('v(2) = w(2) + y(2)');

P_gen = connect(G, Sum1, Sum2, Wp, Wu, {'w(1)', 'w(2)', 'u(1)', 'u(2)'}, {'z1', 'z2', 'v(1)', 'v(2)'});
P_gen = minreal(P_gen);

%% Question 2.7

[K_2, CL_2, gamma_2, info_2] = hinfsyn(P_gen, 2, 2);

L = G * K_2;
gen_nyq_arg = eye(2) + L;
gen_nyq = gen_nyq_arg(1, 1) * gen_nyq_arg(2, 2) - gen_nyq_arg(1, 2) * gen_nyq_arg(2, 1);

figure(3)
nyquist(gen_nyq)

%% Question 2.8

G = tf(FWT);

% Plant
G.InputName = {'u(1)'; 'u(2)'; 'd'};
G.OutputName = 'y';

% Controller
K_2.InputName = 'v';
K_2.OutputName = 'u';

% Sum block
Sum1 = sumblk('v(1) = -w(1) + y(1)');
Sum2 = sumblk('v(2) = -w(2) + y(2)');

P = connect(G, K_2, Sum1, Sum2, {'w(1)', 'w(2)', 'd'}, {'y(1)', 'y(2)'});
P = minreal(P);

figure(4)
step(P)
