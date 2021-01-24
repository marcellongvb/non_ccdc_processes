function [R, S] = ccdc_robustness_outer(W, dims, sdp_prob, robustness, k)
% CCDC ROBUSTNESS MEASURED USING OUTER APPROXIMATION OF CCDC SET
% This function is used for calculating the CCDC robustness of a given
% process W with dimensions (dAI, dAO, dBI).
% Inputs:
% W -- Input process 
% dims -- Vector of dimensions [dAI dAO dBI]
% sdp_prob -- Optimization problem ('primal', 'dual')
% robustness -- options: 'GR', 'WNR'. 'GR' for generalised robustness,
% 'WNR' for white noise robustness.
% intepreter -- ('cvx', 'yalmip')
% k -- level of symmetric extensions
% Outputs:
% R -- Robustness of the process W
% S -- Optimal witness for W (No witness returns from primal problems with
% yalmip)

dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

if nargin < 5
    k = 1;
end

switch sdp_prob
    case 'primal'
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
                R = trace(Omega)/dAO;
                S = full(S);
                
%                 WHITE NOISE ROBUSTNESS PRIMAL
            case 'WNR'
                cvx_begin SDP
                variable R(1,1)
                variable Wdc(dAI*dAO*dBI, dAI*dAO*dBI) complex semidefinite
                variable rho(dAI*dBI, dAI*dBI) complex semidefinite
                dual variable S
                
                minimize R
                
                subject to  
              
                SymmetricExtension(Wdc,k,[dAI dAO*dBI],1,1);
                S: (1-R)*W+R*eye(dAI*dAO*dBI)/(dAI*dBI)==Wdc+PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2],[dAI dBI dAO]);                
                cvx_end
                S = full(S);
        end
    case 'dual'
        switch robustness
            case 'GR'
            yalmip('clear');
            S = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');
            S_dc = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex'); 
            S_dc_AI = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex'); % Component for verifying positive partial transpose on AI

            F = [S_dc + traceandrep(S_dc, [2 3], dims) - traceandrep(S_dc, 3, dims) == 0]; % Orthogonal projection of S_dc_perp onto L_{A->B}
            F = [F, PartialTranspose(S_dc_AI, 1, [dAI dAO*dBI]) >=0]; % Positive Partial Transposition on A_I
            F = [F, PartialTrace(S, 2, dims)>=0]; % CC constraint
            F = [F, S-S_dc_AI+S_dc>=0]; % Witness structure for PPT on A_I and valid ordered process for DC
            F = [F,eye(dAI*dAO*dBI) * (1 + trace(S*W)) - dAO * S >=0]; % Generalized Robustness constraint
            
            R = trace(S*W);        
            ops = sdpsettings('solver', 'mosek', 'verbose', 1);
            SOLUTION=optimize(F,R,ops);        
            R = -value(R);
            S = value(S);               
            case 'WNR'
            yalmip('clear');
            S = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');
            S_dc = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex'); 
            S_dc_AI = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex'); % Component for verifying positive partial transpose on AI

            F = [S_dc + traceandrep(S_dc, [2 3], dims) - traceandrep(S_dc, 3, dims) == 0]; % Orthogonal projection of S_dc_perp onto L_{A->B}
            F = [F, PartialTranspose(S_dc_AI, 1, [dAI dAO*dBI]) >=0]; % Positive Partial Transposition on A_I
            F = [F, PartialTrace(S, 2, dims)>=0]; % CC constraint
            F = [F, S-S_dc_AI+S_dc>=0]; % Witness structure for PPT on A_I and valid ordered process for DC
            F = [F,trace(S)-(dAI*dBI)*trace(S*W)<=(dAI*dBI)]; % White Noise Robustness constraint
            
            R = trace(S*W);        
            ops = sdpsettings('solver', 'mosek', 'verbose', 1);
            SOLUTION=optimize(F,R,ops);        
            R = -value(R);
            S = value(S);     
        end
end
S = (S + S')/2;
end
