function Wopt = FireflyAlgorithm(iterMax, alpha, V_pattern, W0, ~, PdM)
% Firefly Algorithm (FA) following the phase-based procedure described.
%
% Each firefly encodes a phase vector Psi (size |alpha| x 1) with values in
% [-pi, pi]. After optimization, the best Psi is converted to a beamforming
% vector via a least-squares (LS) closed-form, analogous to the last step
% of ILS:
%   p_v = exp(1j * Psi_best)
%   w_opt = (V V^H)^{-1} V D_v p_v, where D_v = diag(Pd)
%
% Interface matches twoStepILS:
%   Wopt = FireflyAlgorithm(iterMax, alpha, V_pattern, W0, PM, PdM)

% --- Setup ---
V = V_pattern(:, alpha);      % M x K
Pd = PdM(:, alpha);
Pd = Pd(:);                   % K x 1
K = length(alpha);

% --- Parameters (as specified) ---
N = 100;       % number of fireflies
beta0 = 1.0;   % base attractiveness
gamma = 1.0;   % light absorption
alphaFA = 0.2; % randomness

% If iterMax is small (e.g., 50 in main.m), FA may underfit vs images.
% Keep iterMax as provided by caller to preserve drop-in compatibility.

% --- Initialize Psi population in [-pi, pi] ---
Psi = -pi + 2*pi*rand(K, N);     % each column is one firefly

% Optional: use W0 to seed one firefly (helps converge)
if ~isempty(W0)
    W0 = W0(:);
    if norm(W0) > 0
        % phase of initial pattern at alpha points
        PdP = (W0' * V) * pinv(diag(Pd));
        Psi(:, 1) = angle(PdP(:));
        Psi(:, 1) = localWrapToPi(Psi(:, 1));
    end
end

costVal = zeros(1, N);
for i = 1:N
    costVal(i) = faCostFromPsi(Psi(:, i), V, Pd);
end

% --- FA iterations (movement) ---
for t = 1:iterMax
    fprintf("FA Iteration %i, Best cost: %f\n", t, min(costVal));
    for i = 1:N
        for j = 1:N
            if costVal(j) < costVal(i)
                rij = norm(Psi(:, i) - Psi(:, j));
                beta = beta0 * exp(-gamma * rij^2);

                Psi(:, i) = Psi(:, i) + beta * (Psi(:, j) - Psi(:, i)) + alphaFA * (rand(K, 1) - 0.5);
                Psi(:, i) = localWrapToPi(Psi(:, i));

                costVal(i) = faCostFromPsi(Psi(:, i), V, Pd);
            end
        end
    end
end

% --- Post-processing (LS) ---
[~, bestIdx] = min(costVal);
psiBest = Psi(:, bestIdx);
pvOpt = exp(1j * psiBest);                 % K x 1
Wopt = (V * V') \ (V * (diag(Pd) * pvOpt)); % M x 1
Wopt = Wopt / max(norm(Wopt), eps);
end

function J = faCostFromPsi(Psi, V, Pd)
% Cost(Psi) based on resulting LS beamformer w(Psi).
pv = exp(1j * Psi);
w = (V * V') \ (V * (diag(Pd) * pv));
p = abs(w' * V).';
err = p - Pd;
J = norm(err, 2)^2;
end

function x = localWrapToPi(x)
% Wrap values to [-pi, pi] without toolboxes.
x = mod(x + pi, 2*pi) - pi;
end
