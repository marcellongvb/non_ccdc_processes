% This code tests the CCDC non-classicality of the fully non-classical
% process (with CNOT decoder)

% Decoder calculation:
CNOT = [1 0 0 0;
        0 1 0 0;
        0 0 0 1;
        0 0 1 0];

Decoder = zeros(16,16);

for i=0:1
    for j=0:1
        for k=0:1
            for l=0:1
                Decoder = Decoder + kron(kron(ket(i),ket(j)) * kron(bra(k),bra(l)),(CNOT * kron(ket(i),ket(j)) * kron(bra(k),bra(l)) * CNOT'));
            end
        end
    end
end

% Decoder Hilbert spaces: Aux Ao/Bi Aux_o
% Removing Aux_o

Decoder = PartialTrace(Decoder,4,[2 2 2 2]);

% Input state:


phiplus = 1/sqrt(2) * (kron(ket(0),ket(0)) + kron(ket(1),ket(1)));
phiplus = phiplus * phiplus';

% phiplus Hilbert spaces: Ai Aux 


% For obtaining the process W, we need to do the link product of phiplus
% with Decoder. I divide this into two steps:
% Step 1: define term_1 as the extension of phiplus to Ao and Bi with
% partial transpose on Aux

term_1 = kron(PartialTranspose(phiplus,2,[2 2]),eye(4));

% term_1 Hilbert spaces: Ai Aux Ao Bi

% Defining term_2 as the extension of Decoder with Ai:

term_2 = kron(eye(2), Decoder);
% term_2 Hilbert spaces: Ai Aux Ao Bi

% Making the link product

W = PartialTrace(term_1 * term_2, 2, [2 2 2 2]);

[J, S, ~] = ccdc_complete([2 2 2 1 1 1], 'GR', 'primal', W)
