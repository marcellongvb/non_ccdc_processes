function [states] = NRandomPureStates(d, N)
% This function returns a set of random pure states with dimension d
states=zeros(d,d,N);
for i=1:N
psi = RandomStateVector(d);
states(:,:,i) = psi*psi';
end

end