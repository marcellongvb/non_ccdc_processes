% This function calculates, for the Process
% Matrix Formalism, any matrix corresponding to partially tracing some
% subsystems and replacing them by normalized identities
% 
% Inputs:
% 
% W: Input process matrix
% d = vector of dimensions (e.g.: [dAI dAO dBI dBO dCI dCO ...])
% part = subsystems to be replaced by normalized identities
% (e.g: AI and BO = [2 4])


function [partW] = traceandrep(W,part,d) 


idxrest = setdiff(1:length(d),part); 
% dimpart = d(part);

dpart = 1;
for i=1:length(part)
dpart = dpart * d(part(i));
i=i+1;
end

[~,neworder] = sort([part idxrest]);

partW = kron2(eye(dpart)/(dpart),PartialTrace(W,part,d));
partW = PermuteSystems(partW,neworder,[d(part) d(idxrest)]);

end
