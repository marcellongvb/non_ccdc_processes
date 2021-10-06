function [Wmax] = process_max_wit_viol(S, dims)
% This function calculates the process Wmax that maximally violates a
% non-classical CCDC witness S
% Inputs:
% S -- Non-classical CCDC witness
% dims -- dimension vector [dAI dAO dBI]
dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

Sbra=S(:)';
cvx_begin SDP
    variable Wmax(dAI*dAO*dBI,dAI*dAO*dBI) complex semidefinite

   % minimize trace(Wmax*S)
    minimize Sbra*Wmax(:)

    subject to
    trace(Wmax)==dAO;
    traceandrep(Wmax,[2 3],dims)==traceandrep(Wmax,3,dims); 
cvx_end       
Wmax = (Wmax + Wmax')/2;
end




