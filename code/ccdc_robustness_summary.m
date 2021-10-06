function ccdc_robustness_summary(W, dims, states, k)
% This function generates all the non-classical CCDC robustnesses of W,
% obtained from the primal and dual problems of all approximating sets.
% Inputs:
%    W - The process which will have the robustnesses measured
%    dims - dimension vector [dAI dAO dBI] or [dAI dAO dBI dBO dCI] for
%    W == W_FB;
%    states - set of random pure states on AI
%    k - level of symmetric extensions


if size(dims,2)==5  
%     Tripartite processes
        [GR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'GR', states);
        [WNR_primal_up, ~] = ccdc_tripartite_robustness_inner(W, dims, 'WNR', states);
        
        [GR_primal_low_sep, ~] = ccdc_tripartite_robustness_outer(W, dims, 'GR', k);
        [WNR_primal_low_sep, ~] = ccdc_tripartite_robustness_outer(W, dims, 'WNR', k);
        
else
%     Bipartite processes
        [GR_primal_up, ~] = ccdc_robustness_inner(W, dims, 'primal', 'GR', states);
        [WNR_primal_up, ~] = ccdc_robustness_inner(W, dims, 'primal', 'WNR', states);
    
        [GR_primal_low_sep, ~] = ccdc_robustness_outer(W, dims,'GR', k);
        [WNR_primal_low_sep, ~] = ccdc_robustness_outer(W, dims,'WNR', k);
  
        [GR_primal_low_ccdc, ~] = ccdc_robustness_outerNEW(W, dims, 'GR', states, k);
        [WNR_primal_low_ccdc, ~] = ccdc_robustness_outerNEW(W, dims, 'WNR', states, k);
end

if size(dims,2)==3
    fprintf('Summary of robustnesses: \n');
    fprintf('Generalized robustness: \n');
    fprintf('Upper-bound : %.4f\n', GR_primal_up);
    fprintf('Lower-bound (CCDC) : %.4f \n', GR_primal_low_ccdc);
    fprintf('Lower-bound (Separability): %.4f \n', GR_primal_low_sep);

    fprintf('White noise robustness: \n');
    fprintf('Upper-bound : %.4f  \n', WNR_primal_up);
    fprintf('Lower-bound (CCDC) : %.4f \n', WNR_primal_low_ccdc);
    fprintf('Lower-bound (Separability): %.4f \n', WNR_primal_low_sep);
elseif size(dims,2)==5

    fprintf('Summary of robustnesses: \n');
    fprintf('Generalized robustness: \n');
    fprintf('Upper-bound: %.4f \n', GR_primal_up);
    fprintf('Lower-bound: %.4f \n', GR_primal_low_sep);

    fprintf('White noise robustness: \n');
    fprintf('Upper-bound: %.4f \n', WNR_primal_up);
    fprintf('Lower-bound: %.4f \n', WNR_primal_low_sep);
end
