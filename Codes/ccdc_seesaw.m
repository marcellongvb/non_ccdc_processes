function [Rmax, Wmax] = ccdc_seesaw(dims, robustness, approximation, random_state_channel, W)
% This function is a heuristic see-saw for finding a bipartite ordered 
% process with high CCDC robustness (either white noise or generalized
% noise).
% Inputs:
% dims -- vector of dimensions [dAI dAO dBI]
% robustness -- The choice of robustness ('GR' or 'WNR')
% approximation -- The choice of approximation for the CCDC set ('inner' or 'outer') 
% random_input (optional, bool) --  If True, the algorithm selects any random 
% ccdc process to start the see-saw. When false, it constructs a process 
% from a random pure state for AI Aux and a random channel from Aux AO to
% BI. Set to false as default.
% W (optional) -- An input process used to start the see-saw method
% Outputs -- 
% 

dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

%We now set the number of random states used for the inner aproximation
if isequal(approximation, 'inner')
    if dBI > 4
        states = sampled_pure_states(dAI, 150);
    else
        states = sampled_pure_states(dAI, 400);
    end    
end


if nargin < 5
    W = random_bipartite_ordered_process(dims, random_state_channel);
end


R = [];
R(1) = 0;
R(2) = 1;
i = 2;
Wmax = W;

while abs(R(i) - R(i-1)) > 0.0001
    i = i+1;  
    switch approximation
        case 'inner'
       [Rmax, Smax] = ccdc_robustness_inner(Wmax, [dAI dAO dBI], 'primal', robustness, states);
        case 'outer'
       [Rmax, Smax] = ccdc_robustness_outer(Wmax, [dAI dAO dBI], robustness, 1);
    end
    R(i) = Rmax;
    Smax = (Smax + Smax')/2;
    [Wmax] = process_max_wit_viol(Smax, [dAI dAO dBI]);
end
disp('number of iterations:')
i-2
Rmax = R(3:size(R,2));

end
