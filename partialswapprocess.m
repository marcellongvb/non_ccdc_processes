% Partial Swap process
% Decoder: 1/sqrt(2) (Id_b|d \otimes Id_f|e + i Id_b|e \otimes ID_f|d)(SWAP)
% in our notation, it is equivalent to
% 1/sqrt(2) ( Id^(Aux Ao /Aux_o Bi) + i Id^(Aux Ao/Bi Aux_o) 
% Construction decoder in Choi representation:

SWAP = [1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1];
clear i
Decoder = 1/sqrt(2) * (eye(4) + i * SWAP)

decoder_choi = zeros(16,16);
for i = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                decoder_choi = decoder_choi + kron( kron(ket(i),ket(j)) * kron(bra(k),bra(l)), Decoder * kron(ket(i),ket(j) * kron(bra(k),bra(l))) * Decoder');
            end
        end
    end
end
clear i

decoder_choi = PartialTrace(decoder_choi,4,[2 2 2 2]);
% decoder_choi Hilbert spaces: Aux Ao/ Bi

% Input state: Ai Aux
phiplus = 1/sqrt(2) * (kron(ket(0),ket(0)) + kron(ket(1), ket(1)));
phiplus = phiplus * phiplus'

term_1 = kron(PartialTranspose(phiplus,2,[2 2]), eye(4)); % Hilbert spaces: Ai Aux Ao Bi
term_2 = kron(eye(2), decoder_choi); % Hilbert spaces: Ai Aux Ao Bi

W = PartialTrace(term_1 * term_2, 2,[2 2 2 2]);

[J,S ] = ccdc_complete([2 2 2 1 1 1], 'GR', 'primal', W);




