clear

cvx_quiet 'True'
non_ccdc_processes_generator;
d=2;


n = 11
[states etaShrinkANAL etaShrinkNUM] = uniform_qubit_states(n);

m = n^2*2
n_states = size(states,3)
etaSHRINK = [etaShrinkANAL etaShrinkNUM]
for i=1:n_states
    states2(:,:,i)=1/etaShrinkNUM*states(:,:,i) + (1-1/etaShrinkNUM)*eye(d)/d;
end
k = 2;
dims=[2 2 2];
sdp_prob='primal';
robustness='GR';

W=W_PPT;
%W=W222;

%[lowerOLD] = ccdc_robustness_outer(W, dims, sdp_prob, robustness, k)
%lowerOLD=lowerOLD

[lowerNEW] = ccdc_robustness_outerNEW(W, dims, sdp_prob, robustness, states2,k);
%lowerOLD____lowerNEW____=[lowerOLD lowerNEW]
lowerNEW = lowerNEW

[upper] = ccdc_robustness_inner(W, dims, sdp_prob, robustness, states);
%lowerOLD____lowerNEW____upper=[lowerOLD lowerNEW upper]
lowerNEW____upper = [lowerNEW upper]