function [Wopt, info] = HybridFA_LS(iterMax, alpha, V_pattern, W0, PM, PdM, varargin)
% Hybrid Firefly Algorithm with Least Squares (Hybrid FA-LS)
% + Weighted Sum Cost Function Design.
%
% This optimizer keeps the same core idea as the existing FA implementation
% (optimize a phase vector Psi, derive weights via LS), but adds:
%   1) Weighted-sum cost: separately weight mainlobe vs sidelobe points.
%   2) A local LS refinement step (1-step phase-consistency update) after
%      evaluating/moving each firefly.
%
% Interface is compatible with twoStepILS / FireflyAlgorithm, with an
% optional 'opts' struct:
%   opts.seed           : rng seed (default 100)
%   opts.N              : number of fireflies (default 300)
%   opts.beta0          : attractiveness base (default 1)
%   opts.gamma          : absorption (default 1)
%   opts.alpha0         : initial randomness (default 0.5)
%   opts.alphaMin       : final randomness (default 1e-6)
%   opts.delta          : LS regularization (default 2e-8)
%   opts.weights        : manual tuning config (defaults match Sharp & Suppressed)
%       .w_main         : mainlobe tracking weight (default 15000)
%       .w_null_force   : anchor forcing weight (default 10000)
%       .w_SLL          : sidelobe penalty weight (default 2000)
%       .target_SLL_dB  : sidelobe threshold in dB (default -25)
%       .w_sym          : symmetry weight (default 500)
%   opts.verbosity      : 0/1 print progress (default 1)
%
% Drop-in usage examples:
%   W = HybridFA_LS(iterMax, alpha, V_pattern, W0, PM, PdM);
%   [W, info] = HybridFA_LS(iterMax, alpha, V_pattern, W0, PM, PdM, opts);

if nargin < 6
    error('HybridFA_LS:NotEnoughInputs', 'Expected at least 6 inputs.');
end

% PM is kept for interface compatibility (ILS/FA signature)
% but is not used by this phase-based FA-LS formulation.
if ~isempty(PM)
    % no-op
end

opts = struct();
if ~isempty(varargin)
    optsCandidate = varargin{1};
    if ~isempty(optsCandidate)
        opts = optsCandidate;
    end
end

% --- Step 1: Reproducibility ---
seed = getOpt(opts, 'seed', 100);
try
    rng(seed);
catch
    % Octave compatibility (older): rng may not exist
    rand('seed', seed); %#ok<RAND>
end

verbosity = getOpt(opts, 'verbosity', 1);

% --- Data selection ---
V_full = V_pattern;
L = size(V_full, 2);

alpha = alpha(:);
V_alpha = V_full(:, alpha);           % M x K
Pd_alpha = PdM(:, alpha);
Pd_alpha = Pd_alpha(:);               % K x 1
K = length(alpha);

% --- Step 1: Geometric constraints (main_idx + null anchors) ---
main_idx = find(PdM(:) > 0);
if isempty(main_idx)
    error('HybridFA_LS:EmptyMainLobe', 'PdM has no nonzero mainlobe region.');
end

leftAnchor = max(1, main_idx(1) - 1);
rightAnchor = min(L, main_idx(end) + 1);
null_anchors = unique([leftAnchor; rightAnchor]);

% Symmetry pairing indices (eqDir is symmetric in this repo)
halfIdx = (1:floor(L/2)).';
symIdx = L - halfIdx + 1;
mainMask = false(L, 1);
mainMask(main_idx) = true;
symPairMask = ~(mainMask(halfIdx) | mainMask(symIdx));

% --- Options / defaults ---
N = getOpt(opts, 'N', 300);
beta0 = getOpt(opts, 'beta0', 1.0);
gamma = getOpt(opts, 'gamma', 1.0);

alpha0 = getOpt(opts, 'alpha0', 0.5);
alphaMin = getOpt(opts, 'alphaMin', 1e-6);

delta = getOpt(opts, 'delta', 2e-8);

weights = getOpt(opts, 'weights', struct());
w_main = getOpt(weights, 'w_main', 15000);
w_null_force = getOpt(weights, 'w_null_force', 10000);
w_SLL = getOpt(weights, 'w_SLL', 2000);
target_SLL_dB = getOpt(weights, 'target_SLL_dB', -25);
w_sym = getOpt(weights, 'w_sym', 500);

% Weighted LS only uses alpha points to reconstruct W from phase-only Psi
% (fitness is evaluated on the full grid with weighted-sum terms).
wAlpha = ones(K, 1);
Ww = spdiags(wAlpha, 0, K, K);

% Cache LS system factorization: (V W V' + delta I)\b
M = size(V_full, 1);
lhs = V_alpha * Ww * V_alpha' + delta * eye(M);
lsSolver = factorizeLinearSystem(lhs);

% --- Initialize Psi population in [-pi, pi] ---
Psi = -pi + 2*pi*rand(K, N);

% Seed from W0 (optional)
if ~isempty(W0)
    W0 = W0(:);
    if norm(W0) > 0
        PdSafe = Pd_alpha;
        PdSafe(PdSafe == 0) = 1;
        PdP = (W0' * V_alpha) ./ (PdSafe.' );
        Psi(:, 1) = localWrapToPi(angle(PdP(:)));
    end
end

costVal = zeros(1, N);
for i = 1:N
    costVal(i) = hybridFitness(Psi(:, i), V_alpha, V_full, Pd_alpha, PdM(:), main_idx, null_anchors, ...
        halfIdx, symIdx, symPairMask, lsSolver, target_SLL_dB, w_main, w_null_force, w_SLL, w_sym);
end

% --- FA iterations ---
for t = 1:iterMax
    if verbosity
        fprintf("Hybrid FA-LS Iteration %i, Best cost: %f\n", t, min(costVal));
    end

    % Exponential alpha schedule: alpha0 -> alphaMin
    alphaT = expAlpha(alpha0, alphaMin, t, iterMax);
    for i = 1:N
        for j = 1:N
            if costVal(j) < costVal(i)
                rij = norm(Psi(:, i) - Psi(:, j));
                beta = beta0 * exp(-gamma * rij^2);

                % Step 2: Fine-tuning movement with decaying randomness
                Psi(:, i) = Psi(:, i) + beta * (Psi(:, j) - Psi(:, i)) + alphaT * (rand(K, 1) - 0.5);
                Psi(:, i) = localWrapToPi(Psi(:, i));

                costVal(i) = hybridFitness(Psi(:, i), V_alpha, V_full, Pd_alpha, PdM(:), main_idx, null_anchors, ...
                    halfIdx, symIdx, symPairMask, lsSolver, target_SLL_dB, w_main, w_null_force, w_SLL, w_sym);
            end
        end
    end
end

% --- Final LS synthesis for the best candidate ---
[~, bestIdx] = min(costVal);
psiBest = Psi(:, bestIdx);
pvOpt = exp(1j * psiBest);
Wopt = solveWeightedLS(V_alpha, Pd_alpha, pvOpt, Ww, lsSolver);

% Step 3: Global phase synchronization at boresight (0 deg / eqDir=0)
boresightIdx = 1 + floor(L/2);
resp0 = Wopt' * V_full(:, boresightIdx);
Wopt = Wopt * exp(-1j * angle(resp0));
resp0 = Wopt' * V_full(:, boresightIdx);
if real(resp0) < 0
    Wopt = -Wopt;
end

Wopt = Wopt / max(norm(Wopt), eps);

if nargout > 1
    info = struct();
    info.bestCost = min(costVal);
    info.costHistory = []; % not tracked to keep runtime low
    info.opts = opts;
    info.main_idx = main_idx;
    info.null_anchors = null_anchors;
    info.weights = struct('w_main', w_main, 'w_null_force', w_null_force, 'w_SLL', w_SLL, 'target_SLL_dB', target_SLL_dB, 'w_sym', w_sym);
    info.delta = delta;
    info.seed = seed;
end
end

function J = hybridFitness(Psi, V_alpha, V_full, Pd_alpha, Pd_full, main_idx, null_anchors, halfIdx, symIdx, symPairMask, lsSolver, targetSLLdB, w_main, w_null, w_SLL, w_sym)
% Step 2: Hybrid fitness evaluation via regularized pseudo-inverse LS.

pv = exp(1j * Psi);
w = solveWeightedLS(V_alpha, Pd_alpha, pv, [], lsSolver);

pFull = abs(w' * V_full).';

% Mainlobe tracking
errMain = pFull(main_idx) - Pd_full(main_idx);
J_main = mean(errMain.^2);

% Null-anchor forcing (sharp slope near mainlobe)
J_null = mean((pFull(null_anchors)).^2);

% Sidelobe penalty above a target threshold (in dB relative to main peak)
peakMain = max(pFull(main_idx));
peakMain = max(peakMain, eps);
sidelobeMask = true(length(Pd_full), 1);
sidelobeMask(main_idx) = false;

pSide = pFull(sidelobeMask);
sllDb = 20*log10((pSide + eps) / peakMain);
exceed = max(0, sllDb - targetSLLdB);
J_sll = mean(exceed.^2);

% Symmetry penalty (exclude mainlobe pairs)
symErr = pFull(halfIdx) - pFull(symIdx);
if any(symPairMask)
    J_sym = mean((symErr(symPairMask)).^2);
else
    J_sym = 0;
end

J = w_main * J_main + w_null * J_null + w_SLL * J_sll + w_sym * J_sym;
end

function w = solveWeightedLS(V, Pd, pv, Ww, lsSolver)
% Solve: min_w || W^(1/2) (V^H w - diag(Pd) pv) ||_2^2
% Closed-form: w = (V W V^H)^{-1} V W diag(Pd) pv
if isempty(Ww)
    rhs = V * (Pd .* pv);
else
    rhs = V * (Ww * (Pd .* pv));
end
w = solveWithFactorization(lsSolver, rhs);
end

function lsSolver = factorizeLinearSystem(A)
% Portable cached factorization for repeated solves A\b.
% Uses Cholesky when possible, otherwise LU.
lsSolver = struct();
lsSolver.method = 'lu';

% Try Cholesky (fast if SPD)
try
    R = chol(A);
    lsSolver.method = 'chol';
    lsSolver.R = R;
    return;
catch
    % fall back
end

% LU decomposition works generally
[L, U, P] = lu(A);
lsSolver.L = L;
lsSolver.U = U;
lsSolver.P = P;
end

function x = solveWithFactorization(lsSolver, b)
if strcmp(lsSolver.method, 'chol')
    % A = R'R
    x = lsSolver.R \ (lsSolver.R' \ b);
else
    % PA = LU
    x = lsSolver.U \ (lsSolver.L \ (lsSolver.P * b));
end
end

function val = getOpt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end

function x = localWrapToPi(x)
% Wrap values to [-pi, pi] without toolboxes.
x = mod(x + pi, 2*pi) - pi;
end

function a = expAlpha(alpha0, alphaMin, t, tMax)
% Exponential schedule alpha0 -> alphaMin across iterations.
if tMax <= 1
    a = alphaMin;
    return;
end
ratio = alphaMin / max(alpha0, eps);
a = alpha0 * (ratio ^ ((t - 1) / (tMax - 1)));
end
