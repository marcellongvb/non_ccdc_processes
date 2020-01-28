% This script generates pseudo-random Process Matrix by sorting a state vector and
% creating a valid positive semidefinite process from it

function [W] = RandomProcessMatrix(d,probtype)

dAI = d(1);
dAO = d(2);
dBI = d(3);
dBO = d(4);
dCI = d(5);
dCO = d(6);

W = RandomStateVector(dAI*dAO*dBI*dBO*dCI*dCO);
W = W * W';

if isequal(probtype,'ccdc')
    W = traceandrep(W,6,d) - traceandrep(W,[5 6],d) + traceandrep(W,[4 5 6],d) - traceandrep(W,[3 4 5 6],d) +  traceandrep(W,[2 3 4 5 6],d);
elseif isequal(probtype,'order')
    if dCI == 1;
       W = traceandrep(W,2,d) + traceandrep(W,4,d) - traceandrep(W,[2 4],d) - traceandrep(W,[3 4],d) + traceandrep(W,[2 3 4],d) ...
            - traceandrep(W,[1 2],d) + traceandrep(W,[1 2 4],d);
    else
       W = validprocessproj(W,d);
    end
end

if eigs(W,1,'sr') <= 0
    W = W - eigs(W,1,'sr')*eye(dAI*dAO*dBI*dBO*dCI*dCO);
    W = W/trace(W) * (dAO*dBO*dCO);
end

end

