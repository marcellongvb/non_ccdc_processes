function braout = bra(i,d)
% Function for returning <i| of dimension d in the computational basis.
% Standard dimension d = 2 (<0|, <1|)
if nargin < 2
    braout = ket(i)';
else
   braout = ket(i,d)';    
end

