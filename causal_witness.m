% This script is a complete tool for witnessing causal nonseparability. 
% It is able to determine whether a bipartite or tripartite process matrix
% is a convex combination of causal processes. Also, if it is a causal
% process, it is able to determine whether it is a convex combination of
% common-cause and direct-cause processes, or if it is something more than
% that.
% INPUTS:
% d: vector of dimensions ([dAI dAO dBI dBO dCI dCO])
% prob: 'order' or 'ccdc', determines when the program will stop entering
% inside the set of valid processes. If the option is 'order', then it
% stops just by checking if a process is a convex combination of AB and BA
% or not. If the option 'ccdc'is set, the program checks if the process
% violates the order witness and, if not, then calculates witness of
% common-cause/direct-cause.
% robustness: 'GR' for Generalized Robustness or 'RR' for Random Robustness
% sdp_prob: 'primal' for witnesses or 'dual' for decomposition of process
% into sum of causal processes.
% W: Process to be tested. If not given, a default process is used
% depending on the chosen scenario

function [varargout] = causal_witness(d,robustness,W)


% Checking if process is causal
    [J_order,~] = order_sdp(d,robustness,'primal',W);
    
    if J_order > -0.0001
        if (J_order) < 0.0001
            disp('The Process W is a Border order process')
        end
        disp('The Process W is causal')
        disp('Checking if it is a Common-cause/direct-cause process')
        
        Labc = traceandrep(W,6,d) - traceandrep(W,[5 6],d) + traceandrep(W,[4 5 6],d) ...
        - traceandrep(W,[3 4 5 6],d) + traceandrep(W,[2 3 4 5 6],d);
        Lbac = traceandrep(W,6,d) - traceandrep(W,[5 6],d) + traceandrep(W,[2 5 6],d) ...
        - traceandrep(W,[1 2 5 6],d) + traceandrep(W,[1 2 4 5 6],d);
        
        if fidelity(Labc,W)/trace(W) > 0.9999
        [J_ccdc,~] = ccdc_complete(d,robustness,'primal',W);
        disp('The process is ordered as A -> B') 
            if J_ccdc > 0.0001
                disp('The process is CCDC Separable')
            elseif J_ccdc < -0.0001
                disp('The process is CCDC Non-separable')
            elseif J_ccdc > -0.0001 && J_ccdc < 0.0001
                disp('The process is a Border CCDC process')
            end
        elseif fidelity(Lbac,W)/trace(W) > 0.9999
             disp('The process is ordered as B -> A') 
            W = PermuteSystems(W,[3 4 1 2 5 6],d); % This reshapes the process into [BI BO AI AO CI CO], so ccdc_complete can work normally
            [J_ccdc,~] = ccdc_complete(d,robustness,'primal',W);
            disp('The process is ordered as A -> B') 
            if J_ccdc > 0.0001
                disp('The process is CCDC Separable')
            elseif J_ccdc < -0.0001
                disp('The process is CCDC Non-separable')
            elseif J_ccdc > -0.0001 && J_ccdc < 0.0001
                disp('The process is a Border CCDC process')
            end
        elseif fidelity(Labc,W)/trace(W) < 0.9999 && fidelity(Lbac,W)/trace(W) < 0.9999
            disp('The Process is a convex combination of A->B and B->A')
        
    elseif J_order < -0.0001
        disp('The Process W is not causal')        
    end
    
    varargout{1} = J_order;
    if J_order > -0.0001
        varargout{2} = J_ccdc;
    else
        varargout{2} = [];
    end        
    
end