% This script generates every robustness for almost every process presented
% in our paper. W222, W224, W_MRSR, W_PPT, W_SEP 

states = zeros(2,2,1000);
for i = 1:1000
    psi = RandomStateVector(2);
    states(:,:,i) = psi*psi';
end

k = 2;

non_ccdc_processes_generator; % This command generates every process mentioned in the paper.

[W222_robustnesses, ~] = ccdc_robustness_summary(W222, [2 2 2], states, k);  

[W224_robustnesses, ~] = ccdc_robustness_summary(W224, [2 2 4], states, k);  

[W_MRSR_robustnesses, ~] = ccdc_robustness_summary(W_MRSR, [2 2 2],states, k);  

[W_PPT_robustnesses, ~] = ccdc_robustness_summary(W_PPT, [2 2 2], states, k);  

[W_SEP_robustnesses, ~] = ccdc_robustness_summary(W_SEP, [2 2 2], states, k);  

fprintf('Robustnesses via primal optimization problem: \n')

fprintf('            |        R_G        |        R_WN        |\n ')
fprintf('W222       |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W222_robustnesses(3), W222_robustnesses(1),  W222_robustnesses(7), W222_robustnesses(5)); 
fprintf('W224       |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W224_robustnesses(3), W224_robustnesses(1),  W224_robustnesses(7), W224_robustnesses(5)); 
fprintf('W_MRSR     |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_MRSR_robustnesses(3), W_MRSR_robustnesses(1),  W_MRSR_robustnesses(7), W_MRSR_robustnesses(5)); 
fprintf('W_PPT      |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_PPT_robustnesses(3), W_PPT_robustnesses(1),  W_PPT_robustnesses(7), W_PPT_robustnesses(5)); 
fprintf('W_SEP      |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_SEP_robustnesses(3), W_SEP_robustnesses(1),  W_SEP_robustnesses(7), W_SEP_robustnesses(5)); 



fprintf('Robustnesses via dual optimization problem: \n')

fprintf('            |        R_G        |        R_WN        \n| ')
fprintf('W222       |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W222_robustnesses(4), W222_robustnesses(2),  W222_robustnesses(8), W222_robustnesses(6)); 
fprintf('W224       |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W224_robustnesses(4), W224_robustnesses(2),  W224_robustnesses(8), W224_robustnesses(6)); 
fprintf('W_MRSR     |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_MRSR_robustnesses(4), W_MRSR_robustnesses(2),  W_MRSR_robustnesses(8), W_MRSR_robustnesses(6)); 
fprintf('W_PPT      |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_PPT_robustnesses(4), W_PPT_robustnesses(2),  W_PPT_robustnesses(8), W_PPT_robustnesses(6)); 
fprintf('W_SEP      |  [%.4f, %.4f] |      [%.4f, %.4f]    | \n ', W_SEP_robustnesses(4), W_SEP_robustnesses(2),  W_SEP_robustnesses(8), W_SEP_robustnesses(6)); 
% end

