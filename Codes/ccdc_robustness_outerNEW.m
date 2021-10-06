function [Rob, S] = ccdc_robustness_outerNEW(W, dims, sdp_prob, robustness, quasi_states, k)
% CCDC ROBUSTNESS MEASURED USING OUTER APPROXIMATION OF CCDC SET
% CONSIDERING TIGHT DIRECT-CAUSE OUTER APPROXIMATION
% This function is used for calculating the CCDC robustness of a given
% process W with dimensions (dAI, dAO, dBI).
% The interpreter used here is CVX.
% Inputs:
% W -- Input process 
% dims -- Vector of dimensions [dAI dAO dBI]
% sdp_prob -- Optimization problem ('primal', 'dual')
% robustness -- options: 'GR', 'WNR'. 'GR' for generalised robustness,
% 'WNR' for white noise robustness.
% quasi_states -- Set of trace-one operators with dimension dAI from which
% the convex hull gives an outer approximation for the set of quantum
% states with dimension dAI
% k -- level of symmetric extensions (redundant constraints of
% entanglement)
% Outputs:
% R -- Robustness of the process W
% S -- Optimal witness for W
% states -- set of random pure states on AI
% Outputs:
% R -- Robustness of the process W
% S -- Optimal witness for W (No witness returns from primal problems with
% yalmip)

dAI = dims(1);
dAO = dims(2);
dBI = dims(3);


N = size(quasi_states,3);


switch sdp_prob
    case 'primal'        
       switch robustness
           
%            PRIMAL GENERALIZED ROBUSTNESS
           case 'GR'
           cvx_begin SDP
           variable Omega(dAI*dAO*dBI, dAI*dAO*dBI) complex semidefinite
           variable Di(dAO*dBI, dAO*dBI, N) complex semidefinite
           variable rho(dAI*dBI, dAI*dBI) complex semidefinite
           dual variable S
          
           
           Wdc = 0;
           for i = 1:N
               Wdc = Wdc + Tensor(quasi_states(:,:,i),Di(:,:,i));
           end
            

           minimize trace(Omega)/dAO
           
           subject to 
           SymmetricExtension(Wdc,k,[dAI dAO*dBI],1,1);
           
    %         Tripartite ordered condition for Omega:          
           PartialTrace(Omega,3,dims)==Tensor(PartialTrace(Omega,[2 3],dims),eye(dAO)/dAO);
    %         Constraint for direct-cause term:
          for i = 1:N
             PartialTrace(Di(:,:,i),2,[dAO dBI])==(eye(dAO)/dAO)*trace(Di(:,:,i));
          end
    %         Generalized robustness contraint:
           S: (1-trace(Omega)/dAO)*W+Omega==Wdc+PermuteSystems(Tensor(rho, eye(dAO)),[1 3 2], [dAI dBI dAO]);
           cvx_end
           Rob = trace(Omega)/(dAO);
           S = full(S);
           
%           PRIMAL WHITE NOISE ROBUSTNESS
           case 'WNR'
           cvx_begin SDP
           variable Di(dAO*dBI, dAO*dBI, N) complex semidefinite
           variable rho(dAI*dBI, dAI*dBI) complex semidefinite
           variable Rob
           dual variable S
           
           Wdc = 0;
           for i = 1:N
               Wdc = Wdc + Tensor(quasi_states(:,:,i),Di(:,:,i));
           end

           minimize Rob
           
           subject to
    %         Constraint for direct-cause term:
          for i = 1:N
             PartialTrace(Di(:,:,i),2,[dAO dBI])==(eye(dAO)/dAO)*trace(Di(:,:,i));
          end
    %         Generalized robustness contraint:
           S: (1-Rob)*W+Rob*(eye(dAI*dAO*dBI)/(dAI*dBI))==Wdc+PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2], [dAI dBI dAO]);
           cvx_end
           S = full(S);
       end
        
    case 'dual'
    switch robustness
        case 'GR'
        yalmip('clear');
            S = sdpvar(dAI*dAO*dBI,dAI*dAO*dBI, 'hermitian', 'complex');
            S_i = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, N, 'hermitian', 'complex');
            
%             Common-cause constraint:
            F = [PartialTrace(S,2,dims)>=0];
            for i=1:N
               F = [F,PartialTrace((S*Tensor(quasi_states(:,:,i),eye(dAI*dBI))+traceandrep(S_i(:,:,i),3,dims)-traceandrep(S_i(:,:,i),[2 3],dims)),1,dims)>=0];
            end
            %             Generalized noise constraint
            F = [F,eye(dAI*dAO*dBI)*(1+trace(S*W))-dAO*S>=0];
            Rob = trace(S*W);
            ops = sdpsettings('solver','mosek', 'verbose', 1);
            SOLUTION = optimize(F,Rob,ops);
            Rob = -value(Rob);
            S = value(S);
        case 'WNR'
            yalmip('clear');
            S = sdpvar(dAI*dAO*dBI,dAI*dAO*dBI, 'hermitian', 'complex');
            S_i = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, N, 'hermitian', 'complex');
            
%             Common-cause constraint:
            F = [PartialTrace(S,2,dims)>=0];
            for i=1:N
               F = [F,PartialTrace((S*Tensor(quasi_states(:,:,i),eye(dAI*dBI))+traceandrep(S_i(:,:,i),3,dims)-traceandrep(S_i(:,:,i),[2 3],dims)),1,dims)>=0];
            end
            %             White noise constraint
            F = [F,trace(S)-(dAI*dBI)*trace(S*W)<=(dAI*dBI)];
            Rob = trace(S*W);
            ops = sdpsettings('solver','mosek', 'verbose', 1);
            SOLUTION = optimize(F,Rob,ops);
            Rob = -value(Rob);
            S = value(S);
        
    end
end
S = (S + S')/2;
end
