function [W_ils, W_fa] = compareILSvsFA(iter_nr, alpha, V_pattern, W0, PM, PdM)
% Compare Iterative Least Squares (ILS) vs Firefly Algorithm (FA)
% This function runs both algorithms and returns their results
% for comparison purposes.
%
% Inputs:
%   iter_nr   - Maximum number of iterations
%   alpha     - Indices of the angles to approximate
%   V_pattern - Array response matrix
%   W0        - Initial weight vector
%   PM        - Initial pattern magnitude
%   PdM       - Desired pattern magnitude
%
% Outputs:
%   W_ils     - Optimized weight vector using ILS
%   W_fa      - Optimized weight vector using FA

fprintf('===== Running Two-Step ILS Algorithm =====\n');
tic;
W_ils = twoStepILS(iter_nr, alpha, V_pattern, W0, PM, PdM);
time_ils = toc;
fprintf('ILS completed in %.4f seconds\n\n', time_ils);

fprintf('===== Running Firefly Algorithm =====\n');
tic;
W_fa = FireflyAlgorithm(iter_nr, alpha, V_pattern, W0, PM, PdM);
time_fa = toc;
fprintf('FA completed in %.4f seconds\n\n', time_fa);

% Compute fitness for comparison
V = V_pattern(:, alpha);

% ILS fitness
P_ils = abs(W_ils' * V);
fitness_ils = mean((P_ils - PdM(:, alpha)).^2);

% FA fitness
P_fa = abs(W_fa' * V);
fitness_fa = mean((P_fa - PdM(:, alpha)).^2);

fprintf('===== Results Comparison =====\n');
fprintf('ILS - Execution time: %.4f s, Final fitness (MSE): %.6f\n', time_ils, fitness_ils);
fprintf('FA  - Execution time: %.4f s, Final fitness (MSE): %.6f\n', time_fa, fitness_fa);

if fitness_ils < fitness_fa
    fprintf('ILS achieved better fitness (lower MSE)\n');
else
    fprintf('FA achieved better fitness (lower MSE)\n');
end

% Plot comparison
eqDir = -1:(1/160):1-(1/160);  % Default equivalent directions
Aq = generateQuantizedArrResponse(size(W0, 1), eqDir);

% Calculate conventional pattern (reference from initial Capon beamforming)
P_refGen = abs(W0' * Aq);

figure('Name', 'ILS vs FA Comparison', 'NumberTitle', 'off');
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
plot(eqDir, 10*log10(P_fa_full/max(P_fa_full)), 'r--', 'LineWidth', 2, 'DisplayName', 'FA');

xlabel('Equivalent directions', 'FontSize', 12);
ylabel('|A|, dB', 'FontSize', 12);
title('Comparison: Two-Step ILS vs Firefly Algorithm', 'FontSize', 14);
xlim([-1 1]);
ylim([-35, 1]);
legend('Location', 'northoutside', 'NumColumns', 4, 'FontSize', 10);
hold off;

% Create a second figure showing the difference
figure('Name', 'Performance Difference', 'NumberTitle', 'off');
subplot(2,1,1);
hold on;
plot(eqDir, abs(P_ils_full/max(P_ils_full) - PdM/max(PdM)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'ILS Error');
plot(eqDir, abs(P_fa_full/max(P_fa_full) - PdM/max(PdM)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'FA Error');
xlabel('Equivalent directions');
ylabel('Absolute Error');
title('Pattern Error Comparison');
legend('Location', 'best');
grid on;
xlim([-1 1]);

subplot(2,1,2);
bar([fitness_ils, fitness_fa]);
set(gca, 'XTickLabel', {'ILS', 'FA'});
ylabel('Mean Squared Error (MSE)');
title('Overall Fitness Comparison');
grid on;

end