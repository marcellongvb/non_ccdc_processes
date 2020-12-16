function [] = bipartite_ordered_test(W,dAI, dAO, dBI)
% This function tests if a given process W satisfies conditions for being a
% valid bipartite ordered process.

% Projecting onto L_{A-> B}
L = W + traceandrep(W, [2 3], [dAI dAO dBI]) - traceandrep(W, 3, [dAI dAO dBI]);

    if min(eig(W)) > -0.00001
        fprintf('The process is positive semidefinite \n');
        if (trace(W) - dAO) <= 0.00001 && trace(W) - dAO > -0.00001
            fprintf('The process satisfies trace(W) = d_AO \n');      
            if max(max(L - W)) < 0.00001 && min(min(L-W)) > -0.00001
                fprintf('The process is in the set of bipartite ordered processes \n \n') 
                fprintf('The process is a valid bipartite ordered process \n')
            else
                fprintf('The process is not in the set L_{A->B} \n')
            end
        else
            fprintf('The operator is not normalized \n')
        end 
    else
    fprintf('The operator W is not positive semidefinite \n')
    end
    
end

