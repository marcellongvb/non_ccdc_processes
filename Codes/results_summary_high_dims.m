% This script generates every robustness for the processes with great 
% dimensions presented in our paper. W_FB, W339

states = zeros(2,2,200);
for i = 1:200
    psi = RandomStateVector(2);
    states(:,:,i) = psi*psi';
end

k = 1;

non_ccdc_processes_generator; % This command generates every process mentioned in the paper.


[GR_FB_upper, ~] = ccdc_tripartite_robustness_inner(W_FB, [2 2 2 2 8], 'GR', states)

[WNR_FB_upper, ~] = ccdc_tripartite_robustness_inner(W_FB, [2 2 2 2 8], 'WNR', states)

[GR_FB_lower, ~] = ccdc_tripartite_robustness_outer(W_FB, [2 2 2 2 8], 'GR', k)

[WNR_FB_lower, ~] = ccdc_tripartite_robustness_outer(W_FB, [2 2 2 2 8], 'WNR', k)

[W_FB_robustnesses, ~] = ccdc_robustness_summary(W_FB, [2 2 2 2 8], states, k);  

[W339_robustnesses, ~] = ccdc_robustness_summary(W339, [3 3 9], states, k);  

fprintf('            |        R_G        |        R_WN        |\n ')
fprintf('W222       |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_FB_robustnesses(2), W_FB_robustnesses(1),  W_FB_robustnesses(4), W_FB_robustnesses(3)); 
fprintf('W224       |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W339_robustnesses(3), W339_robustnesses(1),  W339_robustnesses(7), W339_robustnesses(5)); 
