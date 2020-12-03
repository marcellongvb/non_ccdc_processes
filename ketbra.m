function ketbraout = ketbra(i,j,d)
% Function for returning |i><j| in the computational basis. 
% Standard dimension = 2 (|0><0|, |0><1|, |1><0|, |1><1|)
if nargin < 3
    ketbraout = ket(i) * bra(j);
else
    ketbraout = ket(i,d) * bra(j,d);
end