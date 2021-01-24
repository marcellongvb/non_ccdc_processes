function [robustnesses, order] = ccdc_robustness_summary(W, dims, states, k)
% This function generates all the non-classical CCDC robustnesses of W,
% obtained from the primal and dual problems.
% Inputs:
%    W - The process which will have the robustnesses measured
%    dims - dimension vector [dAI dAO dBI] or [dAI dAO dBI dBO dCI] for
%    W == W_FB;
%    states - set of random pure states on AI
%    k - level of symmetric extensions

if size(dims,2)==5         
    if nargin < 3
        [GR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'GR');
        [WNR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'WNR');
    else
        [GR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'GR', states);
        [WNR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'WNR', states);
    end    
    if nargin < 4
        [GR_primal_low, ~] = ccdc_tripartite_robustness_outer(W, dims, 'GR');
        [WNR_primal_low, ~] = ccdc_tripartite_robustness_outer(W, dims, 'WNR');
    else
        [GR_primal_low, ~] = ccdc_tripartite_robustness_outer(W, dims, 'GR', k);
        [WNR_primal_low, ~] = ccdc_tripartite_robustness_outer(W, dims, 'WNR', k);
    end
else
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
        [GR_primal_low, ~] = ccdc_robustness_outer(W, dims, 'primal', 'GR');
        [WNR_primal_low, ~] = ccdc_robustness_outer(W, dims, 'primal', 'WNR');
    else
        [GR_primal_low, ~] = ccdc_robustness_outer(W, dims, 'primal', 'GR', k);
        [WNR_primal_low, ~] = ccdc_robustness_outer(W, dims, 'primal', 'WNR', k);
    end

    [GR_dual_low, ~] = ccdc_robustness_outer(W, dims, 'dual', 'GR');
    [WNR_dual_low, ~] = ccdc_robustness_outer(W, dims, 'dual', 'WNR');
end

if size(dims,2)==3
    order = ['GR_primal_up, ', 'GR_dual_up, ', 'GR_primal_low, ', 'GR_dual_low, ', 'WNR_primal_up, ', 'WNR_dual_up, ', 'WNR_primal_low, ', 'WNR_dual_low.'];
    robustnesses = [GR_primal_up, GR_dual_up, GR_primal_low, GR_dual_low, WNR_primal_up, WNR_dual_up, WNR_primal_low, WNR_dual_low];

    fprintf('Summary of robustnesses: \n');
    fprintf('Generalized robustness: \n');
    fprintf('Upper-bound - Primal: %.4f, Dual: %.4f. \n', GR_primal_up, GR_dual_up);
    fprintf('Lower-bound - Primal: %.4f, Dual: %.4f. \n', GR_primal_low, GR_dual_low);

    fprintf('White noise robustness: \n');
    fprintf('Upper-bound - Primal: %.4f, Dual: %.4f. \n', WNR_primal_up, WNR_dual_up);
    fprintf('Lower-bound - Primal: %.4f, Dual: %.4f. \n', WNR_primal_low, WNR_dual_low);
elseif size(dims,2)==5
    order = ['GR_primal_up, ', 'GR_primal_low, ', 'WNR_primal_up, ', 'WNR_primal_low, '];
    robustnesses = [GR_primal_up, GR_primal_low, WNR_primal_up, WNR_primal_low];

    fprintf('Summary of robustnesses: \n');
    fprintf('Generalized robustness: \n');
    fprintf('Upper-bound: %.4f \n', GR_primal_up);
    fprintf('Lower-bound: %.4f \n', GR_primal_low);

    fprintf('White noise robustness: \n');
    fprintf('Upper-bound: %.4f \n', WNR_primal_up);
    fprintf('Lower-bound: %.4f \n', WNR_primal_low);
end

% end