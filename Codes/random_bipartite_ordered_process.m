function [W] = random_bipartite_ordered_process(dims, random_state_channel)
% This function generates a random bipartite ordered process W with
% dimensions [dAI dAO dBI].
% Inputs: dims - Vector of dimensions [dAI dAO dBI]
%         random_state_channel (bool): If true, it generates a random
%         process from the link product of a bipartite state and a random
%         channel. If false, it generates a random process from a random
%         density matrix.
dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

switch random_state_channel
    case true        
        
%     State Vector for (AI,aux) (We consider aux and AI having the same
%     dimension)
        d_aux = dAI;
        d_auxO = dAI;

        phi = RandomStateVector(dAI*d_aux);
        rho = phi * phi';
%     Choi Channel from  Aux AO to Bi AuxO
        D = RandomSuperoperator([d_aux*dAO dBI*d_auxO],1,0,0,1);

%     Process [Ai Aux, Ao Bi AuxO] * [Ai, Aux AO BI AuxO]
        term_1 = Tensor(PartialTranspose(rho,2,[dAI d_aux]),eye(dAO*dBI*d_auxO));
        term_2 = Tensor(eye(dAI), D);
        W = PartialTrace(term_1 * term_2, [2 5], [dAI d_aux dAO dBI d_auxO]);
    
    case false
        i = 0;
        W = RandomDensityMatrix(dAI*dAO*dBI);
        W = W + traceandrep(W, [2 3], [dAI dAO dBI]) - traceandrep(W, 3, [dAI dAO dBI]);        
        while min(eig(W)) < 0 && i <= 100
            i = i + 1
    %         Generating random density matrix with dimensions
    %         (dAI*dAO*dBI,dAI*dAO*dBI)
            W = RandomDensityMatrix(dAI*dAO*dBI);
    %         Projecting the density matrix with L_{A->B}(W) = W   
            W = W + traceandrep(W, [2 3], [dAI dAO dBI]) - traceandrep(W, 3, [dAI dAO dBI]);
    %        Normalizing of W
            W = dAO * W;               
        end
         %         Ensuring positive semidefiniteness
        if i ==101 
            lambda = min(eig(W));
            W = dAO * (W - lambda * eye(dAI*dAO*dBI))/trace(W - lambda * eye(dAI*dAO*dBI));
        end        
end
% Ensuring that W is hermitian (removing numerical residuals)
W = (W + W')/2;
end

