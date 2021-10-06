function ketout = ket(i,d)
% Function for returning |i> of dimension d in the computational basis.
% Standard dimension d = 2 (|0>, |1>)

if nargin < 2
    d = 2;
end

ketout = zeros(d,1);
ketout(i+1) = 1;