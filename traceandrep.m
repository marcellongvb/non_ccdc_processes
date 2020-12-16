% This function calculates, for Process Matrices,
% any matrix corresponding to tracing over given parts
% and replacing those parts by normalized identities

% Inputs - W: The matrix in question
%          part = Vector of partitions to be traced
%          d = vector of dimensions 

function [partW] = traceandrep(W, part, d) 

idxrest = setdiff(1:length(d),part); 

dpart = 1;
for i=1:length(part)
dpart = dpart * d(part(i));
end

[~,neworder] = sort([part idxrest]);

partW = Tensor(eye(dpart)/(dpart),PartialTrace(W,part,d));
partW = PermuteSystems(partW,neworder,[d(part) d(idxrest)]);

end
