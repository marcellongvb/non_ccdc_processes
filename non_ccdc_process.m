function [W_out] = non_ccdc_process(process_name)
% This function generates every non-classical CCDC process discussed in our
% work.

% Input- process_name: The name of the process as denoted in our work
%                     (options: 'Wddd2', 'W222', 'W_MRSR', 'W_FB', 'W_SEP',
% 'W_PPT').

    if isequal(process_name,'Wddd2')
        %  Calculating Wddd2 
        d = input('Please inform dimension d=');

        % % % Unnormalized |\phi^+>:
        phiplus = zeros(d^2,1);
        for i=0:(d-1)
             phiplus = phiplus + Tensor(ket(i,d),ket(i,d));
        end

        % % % Channel Ao Bi2:
        Idchoi = phiplus;
        Idchoi = Idchoi * Idchoi';

        % % % State Ai Bi1:
        phiplus = 1/sqrt(d) * phiplus;
        phiplus = phiplus * phiplus';


        % % % Process Ai Bi1 Ao Bi2:
        W_out = Tensor(phiplus,Idchoi);

        % % % Process Ai Ao Bi1 Bi2:
        W_out = PermuteSystems(W_out, [1 3 2 4], [d d d d]);
    
     
    elseif isequal(process_name, 'W222')
    % % % Calculating W222                  
    
    % % % For simplicity, we use the decomposition from Eq.(34) (in terms of
    % Pauli matrices
        W_out = 1/4 * (eye(8) + Tensor(full(Pauli('Z')), eye(2), ...
                full(Pauli('Z'))) + Tensor(full(Pauli('X')), ...
                full(Pauli('X')),full(Pauli('X'))) - ...
                Tensor(full(Pauli('Y')), full(Pauli('X')), ...
                full(Pauli('Y'))));
    elseif isequal(process_name, 'W_MRSR')
    % % % Calculating process W_MRSR
    
        % Decoder: 1/sqrt(2) ( Id^(aux Ao /aux_o Bi) + i SWAP^(aux Ao/aux_o Bi) 
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
        decoder_choi = PartialTrace(decoder_choi,4,[2 2 2 2]);
        % decoder_choi Hilbert spaces: aux Ao/ Bi

        % Input state: Ai aux
        phiplus = 1/sqrt(2) * (Tensor(ket(0),ket(0)) + Tensor(ket(1), ket(1)));
        phiplus = phiplus * phiplus'

        term_1 = Tensor(PartialTranspose(phiplus,2,[2 2]), eye(4)); % Hilbert spaces: Ai aux Ao Bi
        term_2 = Tensor(eye(2), decoder_choi); % Hilbert spaces: Ai aux Ao Bi 

        W_out = PartialTrace(term_1 * term_2, 2,[2 2 2 2]);
    elseif isequal(process_name, 'W_FB')
        % W_FB Choi vector:
        % 1/sqrt(2) ( |0>CI(1) |psi>AIBI |identity>AOCI(2) |identity>BOCI(3) +
        % |1>CI(1)|psi>AICI(2) |identity>AOBI |identity>BOCI(3)) 
        dAI = 2;
        dAO = 2;
        dBI = 2;
        dBO = 2; 
        dCI1 = 2;
        dCI2 = 2;
        dCI3 = 2;
        psi = 1/(sqrt(2)) * (Tensor(ket(0),ket(0)) + Tensor(ket(1),ket(1)));
        phiplus = Tensor(ket(0),ket(0)) + Tensor(ket(1), ket(1));

        term_1 = Tensor(ket(0),psi,phiplus,phiplus);
        term_1 = PermuteSystems(term_1,[2 4 3 6 1 5 7],[dCI1 dAI dBI dAO dCI2 dBO dCI3]);

        term_2 = PermuteSystems(Tensor(ket(1),psi,phiplus,phiplus),[2 4 5 6 1 3 7],[dCI1 dAI dCI2 dAO dBI dBO dCI3]);

        %    W_FB:
        W_out = 1/sqrt(2) * (term_1 + term_2);
        W_out = W_out * W_out';
    elseif isequal(process_name, 'W_SEP')
        plus = (1/sqrt(2)) * (ket(0) + ket(1)); % plus = |+>
        minus = (1/sqrt(2)) * (ket(0) - ket(1)); % minus = |->
        
        W_out = 1/2 * (Tensor(ketbra(0,0),ketbra(0,0) , ketbra(0,0)) + ...
                       Tensor(ketbra(1,1), ketbra(0,0), ketbra(1,1)) + ...                                     
                       Tensor(plus*plus', ketbra(1,1), plus*plus')  + ...
                       Tensor(minus*minus', ketbra(1,1), minus*minus'));
    elseif isequal(process_name, 'W_PPT')
        W_out = 2 * HorodeckiState(0.5, [2 4]);
    end

end
