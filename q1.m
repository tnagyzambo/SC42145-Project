%% Setup

clc;
clear;

load("DATA.mat")

B_siso = B(:, 1);
C_siso = C(1, :);
D_siso = D(1, 1);

sys_siso = ss(A, B_siso, C_siso, D_siso);

G_p = -1 * tf(sys_siso);

%% Question 1.3
s = tf('s');

k_p = 0.15;
tau_i = 0.6;

C_pi = k_p * (tau_i * s + 1) / (tau_i *s);

G_c = C_pi;
sys_open = series(G_c, G_p);
y_r = feedback(sys_open, 1);
c_r = G_c / (1 + G_p * G_c);

bandwidth(y_r)

% Custom bode plot

%Setup
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
w = logspace(-5, 3.4, 1001) * 2 * pi;

[mag_p, phase_p, wout_p] = bode(G_p, w);
[mag_c, phase_c, wout_c] = bode(G_c, w);
[mag_sys, phase_sys, wout_sys] = bode(sys_open, w);

mag_db_p = 20*log10(squeeze(mag_p(1, 1, :)));
phase_p = squeeze(phase_p(1, 1, :));

mag_db_c = 20*log10(squeeze(mag_c(1, 1, :)));
phase_c = squeeze(phase_c(1, 1, :));

mag_db_sys = 20*log10(squeeze(mag_sys(1, 1, :)));
phase_sys = squeeze(phase_sys(1, 1, :));

% Margins
[gm, pm, w_cg, w_cp] = margin(sys_open);
gm = 20 * log10(gm);
w_gc = w_cp;
w_pc = w_cg;

% Time constants
g_tau_i = [1 / tau_i, interp1(w, mag_db_c, 1 / tau_i)];
p_tau_i = [1 / tau_i, interp1(w, phase_c, 1 / tau_i)];

% Plotting
f1 = figure(1);
tiled1 = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'none');

% Gain plot
nexttile(tiled1);
p1 = semilogx(w, mag_db_p, 'color', [0.8500 0.3250 0.0980]); % Plant dynamics
hold on
p2 = semilogx(w, mag_db_c, 'color', [0.9290 0.6940 0.1250]); % Controller dynamics
p3 = semilogx(w, mag_db_sys, 'color', [0 0.4470 0.7410]); % Closed loop system dynamics
yline(0, ':', 'color', [0, 0, 0] + 0.25); % Zero gain line
xlim([w(1), w(end)])

if ~isnan(w_gc)
    xl = xline(w_gc); % Gain crossover frequency
    xl.LineStyle = ':';
    xl.Color = [0, 0, 0] + 0.25;
    xl.Label = sprintf('\\omega_{\\it{gc}} =  %.2fHz', w_gc);
    xl.FontSize = 8;
    xl.LabelVerticalAlignment = 'bottom';
    semilogx(w_gc, 0, '.', 'color', [0, 0, 0] + 0.25);
end

if ~isnan(w_pc)
    xl = xline(w_pc); % Phase crossover frequency
    xl.LineStyle = ':';
    xl.Color = [0, 0, 0] + 0.25;
    xl.Label = sprintf('\\omega_{\\it{pc}} =  %.2fHz', w_pc);
    xl.FontSize = 8;
    xl.LabelVerticalAlignment = 'bottom';
    semilogx(w_pc, 0, '.k');
    semilogx(w_pc, -gm, '.k');
    semilogx([w_pc, w_pc], [0, -gm], 'k');
    text(w_pc, -gm / 2, 'GM ', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontSize', 6)
end

plot(g_tau_i(1), g_tau_i(2), 'k.')
text(g_tau_i(1), g_tau_i(2), '\tau_{I}', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')

hold off
legend([p1(1), p2(1), p3(1)], 'Plant', 'Controller', 'System');
ylabel('Magnitude (dB)');

% Phase plot
nexttile(tiled1);
semilogx(w, phase_p, 'color', [0.8500 0.3250 0.0980]); % Plant dynamics
hold on
semilogx(w, phase_c, 'color', [0.9290 0.6940 0.1250]); % Controller dynamics
semilogx(w, phase_sys, 'color', [0 0.4470 0.7410]); % Closed loop system dynamics
plot(w, 180 * ones(length(w)), ':', 'color', [0, 0, 0] + 0.25); % -180° phase line
ylim([-200, 400])
xlim([w(1), w(end)])

if ~isnan(w_gc)
    xline(w_gc, ':', 'color', [0, 0, 0] + 0.25); % Gain cross over line
    semilogx(w_gc, 180, 'k.');
    semilogx(w_gc, 180 + pm, 'k.');
    text(w_gc, (pm / 2) + 180, ' PM', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 6)
    semilogx([w_gc, w_gc], [180, 180 + pm], 'k');
end

if ~isnan(w_pc)
    xline(w_pc, ':', 'color', [0, 0, 0] + 0.25); % Phase cross over line
    semilogx(w_pc, 180, '.', 'color', [0, 0, 0] + 0.25);
end

plot(p_tau_i(1), p_tau_i(2), 'k.')
text(p_tau_i(1), p_tau_i(2), '\tau_{I}', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
hold off

ylabel('Phase (°)');
xlim([min(w) max(w)])

% Shared
gm_str = sprintf('Gain Margin = \\bf%.1fdB\\rm', gm);
pm_str = sprintf('Phase Margin = \\bf%.1f°\\rm', pm);
[t, s] = title(tiled1, 'Bode Diagram', [gm_str, '   ', pm_str]);
t.FontSize = 16;
s.FontAngle = 'italic';
xlabel(tiled1, 'Frequency (Hz)');
%exportgraphics(gcf, 'images/SC42145_q1_3_bode.png', 'Resolution', 600) 

% Closed-loop step response
% Setup
[step_sys_y, step_sys_t] = step(y_r);
[step_c_y, step_c_t] = step(c_r, step_sys_t);
step_results = stepinfo(step_sys_y, step_sys_t, 'ST', 0.01);
step_st = step_results.SettlingTime;
step_smax = step_results.SettlingMax;
step_sfinal = step_sys_y(end);
step_os = step_results.Overshoot;

% Plot step response
f2 = figure(2);
hsx = plotNy(step_sys_t, {step_sys_y, step_c_y}, [1 2], ...
        'YAxisLabels', {'\omega (rad/s)' '|Cont|'}, ...
        'XAxisLabel', 'Time (s)', ...
        'TitleStr', 'Closed-loop Step Response (From: r To: \omega)', ...
        'Grid', 'on', ...
        'Parent', f2, ...
        'Xlim', [0 step_sys_t(end)], ...
        'Ylim', [0 1.4; 0 6], ...
        'colorord', [0 0.4470 0.7410; 0.9290 0.6940 0.1250], ...
        'lineord', {'-' '-.'}, ...
        'LegendString', {'System' 'Controller'}, ...
        'LegendLoc', 'northeast');

hold on;

% Settling time
if ~isnan(step_st)
    xl = xline(step_st, 'HandleVisibility', 'off');
    xl.LineStyle = ':';
    xl.Color = [0, 0, 0] + 0.25;
    xl.Label = sprintf('T_s =  %.2fs', step_st);
    xl.LabelVerticalAlignment = 'bottom';
end

% Final value & overshoot
if ~isnan(step_sfinal)
    yl = yline(step_sfinal, 'HandleVisibility', 'off');
    yl.LineStyle = ':';
    yl.Color = [0, 0, 0] + 0.25;
    yl.Label = sprintf('y_{final} =  %.2f', step_sfinal);
    yl.LabelVerticalAlignment = 'bottom';
    
    if ~isnan(step_smax) && step_smax > dcgain(y_r)
        step_smax_x = interp1(step_sys_y, step_sys_t, step_smax);
        xline(step_smax_x, ':', 'color', [0, 0, 0] + 0.25, 'HandleVisibility', 'off');
        plot(step_smax_x, dcgain(y_r), '.k', 'HandleVisibility', 'off');
        plot(step_smax_x, step_smax, '.k', 'HandleVisibility', 'off');
        plot([step_smax_x, step_smax_x], [dcgain(y_r), step_smax], 'k', 'HandleVisibility', 'off');
        os_str = sprintf('%%OS = %.2f', step_os);
        text(step_smax_x, step_smax + 0.01, os_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
    end
end

hold off;
%exportgraphics(gcf, 'images/SC42145_q1_3_step.png', 'Resolution', 600)

%% Question 1.4

sys_mimo = FWT;

G_c.u = 'e';
G_c.y = 'Beta (deg)';

Sum = sumblk('e = -r + Omega (rad/s)');

sys_mimo_cl = connect(sys_mimo, G_c, Sum, 'V (m/s)', 'Omega (rad/s)');

% Closed-loop disturbance response
% Setup
[step_sys_d_y, step_sys_d_t] = step(sys_mimo_cl);
step_d_results = stepinfo(step_sys_d_y, step_sys_d_t, 'ST', 0.01);
step_d_st = step_d_results.SettlingTime;
step_d_smax = step_d_results.SettlingMax;
step_d_sfinal = step_sys_d_y(end);
step_d_os = step_d_results.Overshoot;

% Plot step response
f3 = figure(3);
hsx = plotNy(step_sys_d_t, step_sys_d_y, 1, ...
        'YAxisLabels', '\omega (rad/s)', ...
        'XAxisLabel', 'Time (s)', ...
        'TitleStr', 'Closed-loop Disturbance Response (From: V To: \omega)', ...
        'Grid', 'on', ...
        'Parent', f3, ...
        'Xlim', [0 step_sys_d_t(end)+100], ...
        'Ylim', [0 0.5], ...
        'colorord', [0 0.4470 0.7410; 0.9290 0.6940 0.1250], ...
        'lineord', {'-' '-.'}, ...
        'LegendString', {'System'}, ...
        'LegendLoc', 'northeast');

hold on;

% Settling time
if ~isnan(step_d_st)
    xl = xline(step_d_st, 'HandleVisibility', 'off');
    xl.LineStyle = ':';
    xl.Color = [0, 0, 0] + 0.25;
    xl.Label = sprintf('T_s =  %.2fs', step_d_st);
    xl.LabelVerticalAlignment = 'top';
end

% Final value & overshoot
if ~isnan(step_d_sfinal)
    yl = yline(step_d_sfinal, 'HandleVisibility', 'off');
    yl.LineStyle = ':';
    yl.Color = [0, 0, 0] + 0.25;
    yl.Label = sprintf('y_{final} =  %.2f', step_d_sfinal);
    yl.LabelVerticalAlignment = 'bottom';
    
    if ~isnan(step_d_smax) && step_d_smax > dcgain(sys_mimo_cl)
        step_d_smax_x = interp1(step_sys_d_y, step_sys_d_t, step_d_smax);
        xline(step_d_smax_x, ':', 'color', [0, 0, 0] + 0.25, 'HandleVisibility', 'off');
        plot(step_d_smax_x, dcgain(sys_mimo_cl), '.k', 'HandleVisibility', 'off');
        plot(step_d_smax_x, step_d_smax, '.k', 'HandleVisibility', 'off');
        plot([step_d_smax_x, step_d_smax_x], [dcgain(sys_mimo_cl), step_d_smax], 'k', 'HandleVisibility', 'off');
        os_str = sprintf('%%OS = %.2f', step_d_os);
        text(step_d_smax_x, step_d_smax + 0.01, os_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
    end
end

hold off;
%exportgraphics(gcf, 'images/SC42145_q1_4.png', 'Resolution', 600)