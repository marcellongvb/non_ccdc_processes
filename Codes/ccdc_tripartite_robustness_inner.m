function [R, S] = ccdc_tripartite_robustness_inner(W, dims, robustness, states)

% SDP for calculating CCDC robustnesses of tripartite ordered processes 
% with the inner approximation.
% Inputs - W: The process to be measured;
%          dims: dimension vector of the process [dAI dAO dBI dBO dCI] 
%          sdp_prob: type of optimization problem ('primal', 'dual') 
%          robustness: 'GR' for generalized robustness, 'WNR' for white noise
% robustness
%          interpreter: ('cvx', 'yalmip')
%          states: Set of states on A_I (states(dAI,dAI,N)).

% % Process dimensions
    dAI = dims(1);
    dAO = dims(2);
    dBI = dims(3);
    dBO = dims(4);
    dCI = dims(5);

if nargin < 4
    N = 150;
    states=zeros(dAI,dAI,N);
    for i = 1:N
       psi = RandomStateVector(dAI);
       states(:,:,i) = psi * psi';
    end
else
    if size(states,3) == 1
        N = states;
        states = zeros(dAI, dAI, N);
            for i=1:150
                psi = RandomStateVector(dAI);
                states(:,:,i) = psi * psi';
            end
    else
        N = size(states,3);
    end
end
% 
switch robustness
    case 'GR'
        cvx_begin SDP
        variable Omega(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable Wcc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable Wdc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable rho(dAI*dBI, dAI*dBI) complex semidefinite
        variable Di(dAO*dBI, dAO*dBI, N) complex semidefinite
        dual variable S
        
        
        Wdc_sub = 0;
        for i = 1:N
           Wdc_sub = Wdc_sub + Tensor(states(:,:,i),Di(:,:,i)); 
        end
         
        minimize trace(Omega)/(dAO*dBO)
        
        subject to
%         Tripartite ordered conditions for Omega, Wcc and Wdc:
        PartialTrace(Omega, 5, dims)==Tensor(PartialTrace(Omega, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Omega, [3 4 5],dims)==Tensor(PartialTrace(Omega,[2 3 4 5],dims),eye(dAO)/dAO);
        
        PartialTrace(Wcc, 5, dims)==Tensor(PartialTrace(Wcc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wcc, [3 4 5],dims)==Tensor(PartialTrace(Wcc,[2 3 4 5],dims),eye(dAO)/dAO);
        
        PartialTrace(Wdc, 5, dims)==Tensor(PartialTrace(Wdc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wdc, [3 4 5],dims)==Tensor(PartialTrace(Wdc,[2 3 4 5],dims),eye(dAO)/dAO);
        
%         Constraint for N valid channels
        for i=1:N
           PartialTrace(Di(:,:,i),2, [dAO dBI])==trace(Di(:,:,i))*eye(dAO)/dAO; 
        end
        
%         Common-cause constraint:
        PartialTrace(Wcc, [4 5], dims) == PermuteSystems(Tensor(rho, eye(dAO)), [1 3 2], [dAI dBI dAO]);
%         Direct-cause constraint:
        PartialTrace(Wdc, [4 5], dims) == Wdc_sub;
        
%         Generalized robustness constraint
        S: (1-trace(Omega)/(dAO*dBO))*W+Omega==Wcc+Wdc;         
        cvx_end
        R = trace(Omega)/(dAO*dBO);
        S = full(S);
    case 'WNR'
           cvx_begin SDP
        variable Rob
        variable Wcc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable Wdc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable rho(dAI*dBI, dAI*dBI) complex semidefinite
        variable Di(dAO*dBI, dAO*dBI, N) complex semidefinite
        dual variable S
        
        
        Wdc_sub = 0;
        for i = 1:N
           Wdc_sub = Wdc_sub + Tensor(states(:,:,i),Di(:,:,i)); 
        end
         
        minimize Rob
        
        subject to
%         Tripartite ordered conditions for Wcc and Wdc:
        PartialTrace(Wcc, 5, dims)==Tensor(PartialTrace(Wcc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wcc, [3 4 5],dims)==Tensor(PartialTrace(Wcc,[2 3 4 5],dims),eye(dAO)/dAO);
        
        PartialTrace(Wdc, 5, dims)==Tensor(PartialTrace(Wdc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wdc, [3 4 5],dims)==Tensor(PartialTrace(Wdc,[2 3 4 5],dims),eye(dAO)/dAO);
        
%         Constraint for N valid channels
        for i=1:N
           PartialTrace(Di(:,:,i),2, [dAO dBI])==trace(Di(:,:,i))*eye(dAO)/dAO; 
        end
        
%         Common-cause constraint:
        PartialTrace(Wcc, [4 5], dims) == PermuteSystems(Tensor(rho, eye(dAO)), [1 3 2], [dAI dBI dAO]);
%         Direct-cause constraint:
        PartialTrace(Wdc, [4 5], dims) == Wdc_sub;
        
%         Generalized robustness constraint
        S: (1-Rob)*W+Rob*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI)==Wcc+Wdc;         
        cvx_end
        S = full(S);
end
S = (S + S')/2;
end
