function [robustnesses, order] = ccdc_robustness_summary(W, dims, states, k)
% This function generates all the non-classical CCDC robustnesses of W,
% obtained from the primal and dual problems of all approximating sets.
% Inputs:
%    W - The process which will have the robustnesses measured
%    dims - dimension vector [dAI dAO dBI] or [dAI dAO dBI dBO dCI] for
%    W == W_FB;
%    states - set of random pure states on AI
%    k - level of symmetric extensions

dAI = dims(1);

[quasi_states, ~, etaShrinkNUM ] = uniform_qubit_states(11);
N = size(quasi_states,3);
for i=1:N
    quasi_states(:,:,i)= (1/etaShrinkNUM) * quasi_states(:,:,i) + (1-1/etaShrinkNUM) * eye(dAI)/dAI;
end


if size(dims,2)==5  
%     Tripartite processes
    if nargin < 3
        [GR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'GR');
        [WNR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'WNR');
    else
        [GR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'GR', states);
        [WNR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'WNR', states);
    end    
    if nargin < 4
        [GR_primal_low_sep, ~] = ccdc_tripartite_robustness_outer(W, dims, 'GR');
        [WNR_primal_low_sep, ~] = ccdc_tripartite_robustness_outer(W, dims, 'WNR');
    else
        [GR_primal_low_sep, ~] = ccdc_tripartite_robustness_outer(W, dims, 'GR', k);
        [WNR_primal_low_sep, ~] = ccdc_tripartite_robustness_outer(W, dims, 'WNR', k);
    end
else
%     Bipartite processes
    if nargin < 3
        [GR_primal_up, ~] = ccdc_robustness_inner(W, dims, 'primal', 'GR');
        [WNR_primal_up, ~] = ccdc_robustness_inner(W, dims, 'primal', 'WNR');
        [GR_dual_up, ~] = ccdc_robustness_inner(W, dims, 'dual', 'GR');
        [WNR_dual_up, ~] = ccdc_robustness_inner(W, dims, 'dual', 'WNR');
    else
        [GR_primal_up, ~] = ccdc_robustness_inner(W, dims, 'primal', 'GR', states);
        [WNR_primal_up, ~] = ccdc_robustness_inner(W, dims, 'primal', 'WNR', states);
        [GR_dual_up, ~] = ccdc_robustness_inner(W, dims, 'dual', 'GR', states);
        [WNR_dual_up, ~] = ccdc_robustness_inner(W, dims, 'dual', 'WNR', states);
    end

    if nargin < 4
        [GR_primal_low_sep, ~] = ccdc_robustness_outer(W, dims, 'primal', 'GR');
        [WNR_primal_low_sep, ~] = ccdc_robustness_outer(W, dims, 'primal', 'WNR');
    else
        [GR_primal_low_sep, ~] = ccdc_robustness_outer(W, dims, 'primal', 'GR', k);
        [WNR_primal_low_sep, ~] = ccdc_robustness_outer(W, dims, 'primal', 'WNR', k);
    end

    [GR_dual_low_sep, ~] = ccdc_robustness_outer(W, dims, 'dual', 'GR');
    [WNR_dual_low_sep, ~] = ccdc_robustness_outer(W, dims, 'dual', 'WNR');
    
    [GR_primal_low_ccdc, ~] = ccdc_robustness_outerNEW(W, dims, 'primal', 'GR', quasi_states, 2);
    [WNR_primal_low_ccdc, ~] = ccdc_robustness_outerNEW(W, dims, 'primal', 'WNR', quasi_states, 2);
    [GR_dual_low_ccdc, ~] = ccdc_robustness_outerNEW(W, dims, 'dual', 'GR', quasi_states, 2);
    [WNR_dual_low_ccdc, ~] = ccdc_robustness_outerNEW(W, dims, 'dual', 'WNR', quasi_states, 2);
end

if size(dims,2)==3
    order = ['GR_primal_up, ', 'GR_dual_up, ', 'GR_primal_low_sep, ', 'GR_dual_low_sep, ', 'WNR_primal_up, ', 'WNR_dual_up, ', 'WNR_primal_low_sep, ', 'WNR_dual_low_sep, ', 'GR_primal_low_CCDC, ', 'WNR_primal_low_CCDC, ',  'GR_dual_low_CCDC, ', 'WNR_dual_low_CCDC.'];
    robustnesses = [GR_primal_up, GR_dual_up, GR_primal_low_sep, GR_dual_low_sep, WNR_primal_up, WNR_dual_up, WNR_primal_low_sep, WNR_dual_low_sep];
    robustnesses = [robustnesses, GR_primal_low_ccdc, WNR_primal_low_ccdc, GR_dual_low_ccdc, WNR_dual_low_ccdc];

    fprintf('Summary of robustnesses: \n');
    fprintf('Generalized robustness: \n');
    fprintf('Upper-bound - Primal: %.4f, Dual: %.4f. \n', GR_primal_up, GR_dual_up);
    fprintf('Lower-bound (CCDC) - Primal: %.4f, Dual: %.4f. \n', GR_primal_low_ccdc, GR_dual_low_ccdc);
    fprintf('Lower-bound (Separability)- Primal: %.4f, Dual: %.4f. \n', GR_primal_low_sep, GR_dual_low_sep);

    fprintf('White noise robustness: \n');
    fprintf('Upper-bound - Primal: %.4f, Dual: %.4f. \n', WNR_primal_up, WNR_dual_up);
    fprintf('Lower-bound (CCDC) - Primal: %.4f, Dual: %.4f. \n', WNR_primal_low_ccdc, WNR_dual_low_ccdc);
    fprintf('Lower-bound (Separability)- Primal: %.4f, Dual: %.4f. \n', WNR_primal_low_sep, WNR_dual_low_sep);
elseif size(dims,2)==5
    order = ['GR_primal_up, ', 'GR_primal_low, ', 'WNR_primal_up, ', 'WNR_primal_low, '];
    robustnesses = [GR_primal_up, GR_primal_low_sep, WNR_primal_up, WNR_primal_low_sep];

    fprintf('Summary of robustnesses: \n');
    fprintf('Generalized robustness: \n');
    fprintf('Upper-bound: %.4f \n', GR_primal_up);
    fprintf('Lower-bound: %.4f \n', GR_primal_low_sep);

    fprintf('White noise robustness: \n');
    fprintf('Upper-bound: %.4f \n', WNR_primal_up);
    fprintf('Lower-bound: %.4f \n', WNR_primal_low_sep);
end

% end