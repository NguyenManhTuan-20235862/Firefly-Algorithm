function [W_ils, W_fa] = compareILSvsFAandHybrid(iter_nr, alpha, V_pattern, W0, PM, PdM, hybOpts)
% Compare Iterative Least Squares (ILS) vs Firefly Algorithm (FA)
% vs Hybrid Firefly Algorithm with Least Squares (Hybrid FA-LS).
%
% This file intentionally mirrors compareILSvsFA.m so that main.m can switch
% by changing only one line:
%   [W_ref(:,1), ~] = compareILSvsFAandHybrid(...)
%
% Outputs (kept compatible with compareILSvsFA.m):
%   W_ils  - ILS beamformer
%   W_fa   - FA (original) beamformer
%
% Hybrid result is computed and included in plots/printed metrics.

if nargin < 7 || isempty(hybOpts)
    % Defaults follow the "Sharp & Suppressed" strategy described in your
    % Section 5.3 write-up.
    hybOpts = struct();
    hybOpts.seed = 100;
    hybOpts.N = 300;
    hybOpts.beta0 = 1.0;
    hybOpts.gamma = 1.0;
    hybOpts.alpha0 = 0.5;
    hybOpts.alphaMin = 1e-6;
    hybOpts.delta = 2e-8;
    hybOpts.verbosity = 1;
    hybOpts.weights = struct(...
        'w_main', 15000, ...
        'w_null_force', 10000, ...
        'w_SLL', 2000, ...
        'target_SLL_dB', -25, ...
        'w_sym', 500);
end

fprintf('===== Running Two-Step ILS Algorithm =====\n');
tic;
W_ils = twoStepILS(iter_nr, alpha, V_pattern, W0, PM, PdM);
time_ils = toc;
fprintf('ILS completed in %.4f seconds\n\n', time_ils);

fprintf('===== Running Firefly Algorithm (Original) =====\n');
tic;
W_fa = FireflyAlgorithm(iter_nr, alpha, V_pattern, W0, PM, PdM);
time_fa = toc;
fprintf('FA completed in %.4f seconds\n\n', time_fa);

fprintf('===== Running Hybrid FA-LS =====\n');
tic;
W_hyb = HybridFA_LS(iter_nr, alpha, V_pattern, W0, PM, PdM, hybOpts);
time_hyb = toc;
fprintf('Hybrid FA-LS completed in %.4f seconds\n\n', time_hyb);

% Compute fitness for comparison (same as compareILSvsFA.m)
V = V_pattern(:, alpha);

P_ils = abs(W_ils' * V);
fitness_ils = mean((P_ils - PdM(:, alpha)).^2);

P_fa = abs(W_fa' * V);
fitness_fa = mean((P_fa - PdM(:, alpha)).^2);

P_hyb = abs(W_hyb' * V);
fitness_hyb = mean((P_hyb - PdM(:, alpha)).^2);

fprintf('===== Results Comparison =====\n');
fprintf('ILS        - Execution time: %.4f s, Final fitness (MSE): %.6f\n', time_ils, fitness_ils);
fprintf('FA (orig)  - Execution time: %.4f s, Final fitness (MSE): %.6f\n', time_fa, fitness_fa);
fprintf('Hybrid FA-LS - Execution time: %.4f s, Final fitness (MSE): %.6f\n', time_hyb, fitness_hyb);

% Plot comparison (kept consistent with compareILSvsFA.m)
eqDir = -1:(1/160):1-(1/160);  % Default equivalent directions
Aq = generateQuantizedArrResponse(size(W0, 1), eqDir);

% Calculate conventional pattern (reference from initial Capon beamforming)
P_refGen = abs(W0' * Aq);

figure('Name', 'ILS vs FA vs Hybrid Comparison', 'NumberTitle', 'off');
hold on;
grid on;

% Plot desired pattern
plot(eqDir, 10*log10(PdM/max(PdM)), 'm-*', 'LineWidth', 1.5, 'DisplayName', 'Desired');

% Plot conventional 12-element ULA (reference)
plot(eqDir, 10*log10(P_refGen/max(P_refGen)), '--black', 'LineWidth', 1.5, 'DisplayName', 'Conventional 12-element ULA');

% Plot ILS result
P_ils_full = abs(W_ils' * Aq);
plot(eqDir, 10*log10(P_ils_full/max(P_ils_full)), 'b-', 'LineWidth', 2, 'DisplayName', 'ILS');

% Plot FA result
P_fa_full = abs(W_fa' * Aq);
plot(eqDir, 10*log10(P_fa_full/max(P_fa_full)), 'r--', 'LineWidth', 2, 'DisplayName', 'FA (orig)');

% Plot Hybrid result
P_hyb_full = abs(W_hyb' * Aq);
plot(eqDir, 10*log10(P_hyb_full/max(P_hyb_full)), 'g-.', 'LineWidth', 2, 'DisplayName', 'Hybrid FA-LS');

xlabel('Equivalent directions', 'FontSize', 12);
ylabel('|A|, dB', 'FontSize', 12);
title('Comparison: Two-Step ILS vs FA (orig) vs Hybrid FA-LS', 'FontSize', 14);
xlim([-1 1]);
ylim([-35, 1]);
legend('Location', 'northoutside', 'NumColumns', 4, 'FontSize', 10);
hold off;

% Create a second figure showing the difference (add Hybrid)
figure('Name', 'Performance Difference (ILS vs FA vs Hybrid)', 'NumberTitle', 'off');
subplot(2,1,1);
hold on;
plot(eqDir, abs(P_ils_full/max(P_ils_full) - PdM/max(PdM)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'ILS Error');
plot(eqDir, abs(P_fa_full/max(P_fa_full) - PdM/max(PdM)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'FA Error');
plot(eqDir, abs(P_hyb_full/max(P_hyb_full) - PdM/max(PdM)), 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Hybrid Error');
xlabel('Equivalent directions');
ylabel('Absolute Error');
title('Pattern Error Comparison');
legend('Location', 'best');
grid on;
xlim([-1 1]);

subplot(2,1,2);
bar([fitness_ils, fitness_fa, fitness_hyb]);
set(gca, 'XTickLabel', {'ILS', 'FA', 'Hybrid'});
ylabel('Mean Squared Error (MSE)');
title('Overall Fitness Comparison');
grid on;

end
