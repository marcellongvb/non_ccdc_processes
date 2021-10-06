function [Rob, S] = ccdc_robustness_outerNEW(W, dims, robustness, states, k)
% BE CAREFUL This function only works if Alice's input space has dimension 2
% CCDC ROBUSTNESS MEASURED USING OUTER APPROXIMATION OF CCDC SET
% CONSIDERING TIGHT DIRECT-CAUSE OUTER APPROXIMATION
% This function is used for calculating the CCDC robustness of a given
% process W with dimensions (dAI, dAO, dBI).
% The interpreter used here is CVX.
% Inputs:
% W -- Input process 
% dims -- Vector of dimensions [dAI dAO dBI]
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
% S -- Optimal witness for W 
dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

N = size(states,3); %Number of states

%The first step here is to evaluate the quasi-states
%First, we evaluate the srhinking factor, and then we inflate the states
%such that the quasi-states form an outer aproximation
for i=1:N
    vectors(:,i)=Rho2Bloch(states(:,:,i));
end
    [A,b,Aeq,beq]=vert2lcon(vectors');
    etaShrinkNUM=max(b);
quasi_states=NaN(2,2,N);
for i=1:N
    quasi_states(:,:,i)= (1/etaShrinkNUM) * states(:,:,i) + (1-1/etaShrinkNUM) * eye(dAI)/dAI;
end

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
               SymmetricExtension(Wdc,k,[dAI dAO*dBI],1,1);
          for i = 1:N
             PartialTrace(Di(:,:,i),2,[dAO dBI])==(eye(dAO)/dAO)*trace(Di(:,:,i));
          end
    %         Generalized robustness contraint:
           S: (1-Rob)*W+Rob*(eye(dAI*dAO*dBI)/(dAI*dBI))==Wdc+PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2], [dAI dBI dAO]);
           cvx_end
           S = full(S);
       end
        
S = (S + S')/2;
end
