% This code generates the process matrix that maximally violates a witness
% of causal nonseparability for the order scenario.

function [W,J] = causalwitmaxviol(d,S)

dAI = d(1);
dAO = d(2);
dBI = d(3);
dBO = d(4);
dCI = d(5);
dCO = d(6);

yalmip('clear')

if dCI == 1
%   Bipartite Scenario
    W = sdpvar(dAI*dAO*dBI*dBO,dAI*dAO*dBI*dBO,'hermitian','complex');
    Lv = traceandrep(W,2,d) + traceandrep(W,4,d) - traceandrep(W,[2 4],d) - traceandrep(W,[3 4],d) + traceandrep(W,[2 3 4],d) ...
            - traceandrep(W,[1 2],d) + traceandrep(W,[1 2 4],d);
    
    F = [W>=0,W-Lv==0,trace(W)-dAO*dBO==0];
    
    J=real(trace(S*W));

    ops = sdpsettings('solver','mosek','verbose',1);
    ops.mosek.MSK_IPAR_NUM_THREADS=6;
    SOLUTION=optimize(F,J,ops)
    
    W = value(W);
    J = value(J);

   
else
%   Tripartite Scenario
    W = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO,'hermitian','complex');
    Lv = validprocessproj(W,d);
    
    F = [W>=0,W-Lv==0,trace(W)-dAO*dBO*dCO==0];
    
    J=real(trace(S*W));

    ops = sdpsettings('solver','mosek','verbose',1);
    ops.mosek.MSK_IPAR_NUM_THREADS=6;
    SOLUTION=optimize(F,J,ops)
    
    W = value(W);
    J = value(J);
end



end