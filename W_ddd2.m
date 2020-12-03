function [W_ddd2] = W_ddd2(d)

% This function generates the process W_{ddd^2} proposed in our work
% Optional input - d: basis dimension for subsystems 
if nargin < 1
    d = 2;
end

phiplus = zeros(d^2,1);
for i=0:(d-1)
     phiplus = phiplus + kron(ket(i,d),ket(i,d));
end

% % % Channel Ao Bi2
Idchoi = phiplus;
Idchoi = Idchoi * Idchoi';

% % % State Ai Bi1
phiplus = 1/sqrt(d) * phiplus;
phiplus = phiplus * phiplus';


% % Process Ai Bi1 Ao Bi2
W_ddd2 = kron(phiplus,Idchoi);

% % % Process Ai Ao Bi1 Bi2
W_ddd2 = PermuteSystems(W_ddd2, [1 3 2 4], [d d d d]);
end
