function [Rob, S] = ccdc_robustness_outer(W, dims, robustness, k)
% CCDC ROBUSTNESS MEASURED USING OUTER APPROXIMATION OF CCDC SET
% This function is used for calculating the CCDC robustness of a given
% process W with dimensions (dAI, dAO, dBI).
% Inputs:
% W -- Input process 
% dims -- Vector of dimensions [dAI dAO dBI]
% robustness -- options: 'GR', 'WNR'. 'GR' for generalised robustness,
% 'WNR' for white noise robustness.
% intepreter -- ('cvx', 'yalmip')
% k -- level of symmetric extensions
% Outputs:
% R -- Robustness of the process W
% S -- Optimal witness for W 

dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

        switch robustness
%             GENERALIZED ROBUSTNESS PRIMAL
            case 'GR'
                cvx_begin SDP
                variable Omega(dAI*dAO*dBI, dAI*dAO*dBI) complex semidefinite
                variable Wdc(dAI*dAO*dBI, dAI*dAO*dBI) complex semidefinite
                variable rho(dAI*dBI, dAI*dBI) complex semidefinite
                dual variable S
                
                minimize trace(Omega)/dAO
                
                subject to
                
              
                SymmetricExtension(Wdc,k,[dAI dAO*dBI],1,1);                
                PartialTrace(Omega,3,dims)==Tensor(PartialTrace(Omega,[2 3],dims),eye(dAO)/dAO);
                S: (1-trace(Omega)/dAO)*W+Omega==Wdc+PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2],[dAI dBI dAO]);                
                cvx_end
                Rob = trace(Omega)/dAO;
                S = full(S);
                
%                 WHITE NOISE ROBUSTNESS PRIMAL
            case 'WNR'
                cvx_begin SDP
                variable Rob
                variable Wdc(dAI*dAO*dBI, dAI*dAO*dBI) complex semidefinite
                variable rho(dAI*dBI, dAI*dBI) complex semidefinite
                dual variable S
                
                minimize Rob
                
                subject to  
              
                SymmetricExtension(Wdc,k,[dAI dAO*dBI],1,1);
                S: (1-Rob)*W+Rob*eye(dAI*dAO*dBI)/(dAI*dBI)==Wdc+PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2],[dAI dBI dAO]);                
                cvx_end
                S = full(S);
        end

S = (S + S')/2;
end
