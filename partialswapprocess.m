function [W_mrsr] = partialswapprocess()
% This function generates the process W_{MRSR} presented in 
% Ref. MacLean, J. P. W., Ried, K., Spekkens, R. W., & Resch, K. J. (2017).
% Quantum-coherent mixtures of causal relations.
% Nature communications, 8(1), 1-10.



% Decoder: 1/sqrt(2) (Id_b|d \otimes Id_f|e + i Id_b|e \otimes ID_f|d)(SWAP)
% In our notation, it is equivalent to
% 1/sqrt(2) ( Id^(aux Ao /aux_o Bi) + i SWAP^(aux Ao/aux_o Bi) 
% Construction decoder in Choi representation:

SWAP = [1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1];
Decoder = 1/sqrt(2) * (eye(4) + 1i * SWAP);

decoder_choi = zeros(16,16);
for h = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                decoder_choi = decoder_choi + Tensor( Tensor(ket(h),ket(j)) * Tensor(bra(k),bra(l)), Decoder * Tensor(ket(h),ket(j) * Tensor(bra(k),bra(l))) * Decoder');
            end
        end
    end
end
% Change 3 to 4 to test after finishing the code
decoder_choi = PartialTrace(decoder_choi,4,[2 2 2 2]);
% decoder_choi Hilbert spaces: aux Ao/ Bi

% Input state: Ai aux
phiplus = 1/sqrt(2) * (Tensor(ket(0),ket(0)) + Tensor(ket(1), ket(1)));
phiplus = phiplus * phiplus'

term_1 = Tensor(PartialTranspose(phiplus,2,[2 2]), eye(4)); % Hilbert spaces: Ai aux Ao Bi
term_2 = Tensor(eye(2), decoder_choi); % Hilbert spaces: Ai aux Ao Bi 

W_mrsr = PartialTrace(term_1 * term_2, 2,[2 2 2 2]);
end


