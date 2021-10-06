function [R, S] = ccdc_tripartite_robustness_outer(W, dims, robustness, k)

% SDP for calculating CCDC robustnesses of the process W_FB using the inner
% or outer approximation.
% Input - robustness: 'GR' for generalized robustness, 'WNR' for white noise
% robustness
%       - approximation: 'inner' or 'outer'

% % Process dimensions
dAI = dims(1);
dAO = dims(2);
dBI = dims(3);
dBO = dims(4);
dCI = dims(5);

if nargin < 4
    k = 1;
end

switch robustness
%     PRIMAL GENERALIZED ROBUSTNESS:
    case 'GR'
        cvx_begin SDP
        variable Omega(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable Wdc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable Wcc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable rho(dAI*dBI, dAI*dBI)  complex semidefinite
        dual variable S

        minimize trace(Omega)/(dAO*dBO)

        subject to 
        
%         Tripartite ordered conditions for Omega, Wcc and Wdc:
        PartialTrace(Omega, 5, dims)==Tensor(PartialTrace(Omega, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Omega, [3 4 5],dims)==Tensor(PartialTrace(Omega,[2 3 4 5],dims),eye(dAO)/dAO);
        
        PartialTrace(Wcc, 5, dims)==Tensor(PartialTrace(Wcc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wcc, [3 4 5],dims)==Tensor(PartialTrace(Wcc,[2 3 4 5],dims),eye(dAO)/dAO);
        
        PartialTrace(Wdc, 5, dims)==Tensor(PartialTrace(Wdc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wdc, [3 4 5],dims)==Tensor(PartialTrace(Wdc,[2 3 4 5],dims),eye(dAO)/dAO);
        
%         Constraint for common-cause term:
        PartialTrace(Wcc,[4 5],dims)==PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2],[dAI dBI dAO]);
%         Constraint for direct-cause term:
       
        SymmetricExtension(PartialTrace(Wdc,[4 5],dims),k,[dAI dAO*dBI],1,1); 
%         Generalized robustness contraint:
        S: (1-trace(Omega)/(dAO*dBO))*W+Omega==Wdc+Wcc;        
        cvx_end
        R = trace(Omega)/(dAO*dBO);
    
%     PRIMAL WHITE NOISE ROBUSTNESS        
    case 'WNR'
        cvx_begin SDP
        variable Wdc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable Wcc(dAI*dAO*dBI*dBO*dCI, dAI*dAO*dBI*dBO*dCI) complex semidefinite
        variable rho(dAI*dBI, dAI*dBI)  complex semidefinite
        variable R(1,1)
        dual variable S

        minimize R

        subject to 
        
%         Tripartite ordered conditions for Wcc and Wdc:        
        PartialTrace(Wcc, 5, dims)==Tensor(PartialTrace(Wcc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wcc, [3 4 5],dims)==Tensor(PartialTrace(Wcc,[2 3 4 5],dims),eye(dAO)/dAO);
        
        PartialTrace(Wdc, 5, dims)==Tensor(PartialTrace(Wdc, [4 5], dims),eye(dBO)/dBO);
        PartialTrace(Wdc, [3 4 5],dims)==Tensor(PartialTrace(Wdc,[2 3 4 5],dims),eye(dAO)/dAO);
        
%         Constraint for common-cause term:
        PartialTrace(Wcc, [4 5], dims)==PermuteSystems(Tensor(rho, eye(dAO)),[1 3 2], [dAI dBI dAO]);
%         Constraint for direct-cause term:
        SymmetricExtension(PartialTrace(Wdc,[4 5],dims),k,[dAI dAO*dBI],1,1);
%         Generalized robustness contraint:
        S: (1-R)*W+R*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI)==Wdc+Wcc;
        
        cvx_end
end
S = (S + S')/2;
S = full(S);
end
