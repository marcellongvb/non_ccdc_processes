% SDP for testing if a given process W is a convex combination of a
% direct-cause process and a common-cause process.
% input d: vector of dimensions
% 
% input robustness: 'GR' for generalized robustness, 'RR' for random
% robustness
% 
% input sdp_prob: 'primal' (for witness) or 'dual' (for decomposition)
% 
% opt input W: Process matrix to be tested. default: Wcoh

function [varargout] = ccdc_complete_robust(d,robustness,W,CS1,CS2,NA)

    % Input dimensions
    dAI = d(1);
    dAO = d(2);
    dBI = d(3);
    dBO = d(4);
    dCI = d(5);
    dCO = d(6);
% tic
    if nargin < 3
        % Coherent ccdc state vector:
        % 1/sqrt(2) ( |0>CI(1) |psi>AIBI |identity>AOCI(2) |identity>BOCI(3) +
        % |1>CI(1)|psi>AICI(2) |identity>AOBI |identity>BOCI(3)) 

        zero = [1; 0];
        one = [0; 1];

        dCI1 = 2;
        dCI2 = 2;
        dCI3 = 2;
        psi = 1/(sqrt(2)) * (kron(zero,zero) + kron(one,one));
        phiplus = kron(zero,zero) + kron(one, one);

        term_1 = kron2(zero,psi,phiplus,phiplus);
        term_1 = PermuteSystems(term_1,[2 4 3 6 1 5 7],[dCI1 dAI dBI dAO dCI2 dBO dCI3]);

        term_2 = PermuteSystems(kron2(one,psi,phiplus,phiplus),[2 4 5 6 1 3 7],[dCI1 dAI dCI2 dAO dBI dBO dCI3]);

        %    W Switch:
        W = 1/sqrt(2) * (term_1 + term_2);
        W = W * W';
        % W = term_1 * term_1';
        clear term_1 term_2 psi phiplus zero one
    end
    
    if nargin < 4
        CS1=[];
        CS2=[];
    end

    yalmip('clear');

%     if isequal(sdp_prob,'primal')
        S = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
        Sp2 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
        Sp3 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
        S_r = sdpvar(dAI*dAO*dBI,dAI*dAO*dBI, 'hermitian', 'complex');
        S_ort2 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
        S_ort3 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
        S_ort4 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');


        LabcS2 = S_ort2 - traceandrep(S_ort2,5,d) + traceandrep(S_ort2,[4 5],d) - traceandrep(S_ort2,[3 4 5],d) +  traceandrep(S_ort2,[2 3 4 5],d);
        LabcS3 = S_ort3 - traceandrep(S_ort3,5,d) + traceandrep(S_ort3,[4 5],d) - traceandrep(S_ort3,[3 4 5],d) +  traceandrep(S_ort3,[2 3 4 5],d);

        LccS4 = S_ort4 - traceandrep(S_ort4,5,d) + traceandrep(S_ort4,[4 5],d) - traceandrep(S_ort4,[3 4 5],d) +  traceandrep(S_ort4,[2 3 4 5],d) ...
            - traceandrep(S_ort4,[4 5],d) + traceandrep(S_ort4,[3 4 5],d) -  traceandrep(S_ort4,[2 3 4 5],d) + traceandrep(S_ort4,[2 4 5],d);
        % LdcQ2 = Q2 - traceandrep(Q2,5,d) + traceandrep(Q2,[4 5],d) - traceandrep(Q2,[3 4 5],d) +  traceandrep(Q2,[2 3 4 5],d);

        F = [Sp2>=0,Sp3>=0];


        F = [F,LabcS2==0,LabcS3==0,LccS4==0];
        
%         CHECK CONVERGENCE (LAST TESTED: 10.000 AND 50.000
%         violation for 30.000: -0,8232, for 40.000: -0,8231, for 50.000:
%         -0.8230
%         NA  = 1000;
        Wit = kron(S_r,eye(dBO*dCI*dCO)/(dBO*dCI*dCO));
        % dim = [dAI dAO*dBI];
        
        % ==============================
        % canonical constraints - begin
        IdAOdBI = eye(dAO*dBI);
        for jNA=1:(dAO*dBI)
        %   Direct-cause:
            ket1 = kron2(eye(dAI),IdAOdBI(:,jNA),eye(dBO*dCI));
            F = [F,ket1'*Wit*ket1>=0];  
        end
        IdAI = eye(dAI);
        for jNA=1:dAI
        %   Direct-cause:
            ket2 = kron2(IdAI(:,jNA),eye(dAO*dBI*dBO*dCI));
            F = [F,dBI*ket2'*Wit*ket2>=0];  
        end
        % canonical constraints - end
        % ==============================
        
        % ==============================
        % sifted constraints - begin
        for jNA=1:length(CS1)
        %   Direct-cause:
            ket1 = kron2(eye(dAI),CS1{jNA},eye(dBO*dCI));
            F = [F,ket1'*Wit*ket1>=0];  
        end
        for jNA=1:length(CS2)
        %   Direct-cause:
            ket2 = kron2(CS2{jNA},eye(dAO*dBI*dBO*dCI));
            F = [F,dBI*ket2'*Wit*ket2>=0];  
        end
        % sifted constraints - end
        % ==============================

        for jNA=1:NA
        %   Direct-cause:
            C = RandomSuperoperator([dAO dBI],1,0,0,1);
            [ket1,~] = eigs(C,1,'lm');
            ket1 = kron2(eye(dAI),ket1,eye(dBO*dCI));
            ket2 = kron2(RandomStateVector(dAI),eye(dAO*dBI*dBO*dCI));
            F = [F,ket1'*Wit*ket1>=0,dBI*ket2'*Wit*ket2>=0];  
        end
        F = [F,Wit+Sp2+S_ort2-S==0,Sp3+S_ort4+S_ort3-S==0];


        if isequal(robustness,'GR')
            Sp1  = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
            S_ort1 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');

            LabcS1 = S_ort1 - traceandrep(S_ort1,5,d) + traceandrep(S_ort1,[4 5],d) - traceandrep(S_ort1,[3 4 5],d) +  traceandrep(S_ort1,[2 3 4 5],d);

            F = [F,Sp1>=0,LabcS1==0,(eye(dAI*dAO*dBI*dBO*dCI)/(dAO*dBO))-S-Sp1-S_ort1==0]; % Generalized Robustness
        elseif isequal(robustness,'RR')
            F = [F,trace(S*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI))==1]; % Random Robustness
        end
        J=real(trace(S*W));

        ops = sdpsettings('solver','mosek','verbose',1);
        ops.mosek.MSK_IPAR_NUM_THREADS=6;
        SOLUTION=optimize(F,J,ops)
%         t = toc;
%         disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))

        varargout{1} = value(J);
        varargout{2} = value(S);
        varargout{3} = value(Wit);
%         varargout{3} = W;
        
%     PENSAR EM COMO SERIA O DUAL NESSE CASO
%     elseif isequal(sdp_prob,'dual')
%         Wcc = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
%         Wdc = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
% 
%         Wdc_transp = PartialTranspose(traceandrep(Wdc,[4 5],d),1,d);
%         LabcWdc = Wdc - traceandrep(Wdc,5,d) + traceandrep(Wdc,[4 5],d) - traceandrep(Wdc,[3 4 5],d) +  traceandrep(Wdc,[2 3 4 5],d);
%         LabcWcc = Wcc - traceandrep(Wcc,5,d) + traceandrep(Wcc,[4 5],d) - traceandrep(Wcc,[3 4 5],d) +  traceandrep(Wcc,[2 3 4 5],d);
% 
% 
%         LccWcc = Wcc - traceandrep(Wcc,5,d) + traceandrep(Wcc,[4 5],d) - traceandrep(Wcc,[3 4 5],d) +  traceandrep(Wcc,[2 3 4 5],d) ...
%         - traceandrep(Wcc,[4 5],d) + traceandrep(Wcc,[3 4 5],d) -  traceandrep(Wcc,[2 3 4 5],d) + traceandrep(Wcc,[2 4 5],d);
% 
% 
%         F = [Wdc>=0,Wcc>=0,Wdc_transp>=0,Wcc-LabcWcc==0,Wdc-LabcWdc==0,Wcc-LccWcc==0];
% 
%         if isequal(robustness,'GR')
%             Omega = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
% 
%             LabcOmega = Omega - traceandrep(Omega,5,d) + traceandrep(Omega,[4 5],d) - traceandrep(Omega,[3 4 5],d) +  traceandrep(Omega,[2 3 4 5],d);
% 
%             F = [F,Omega>=0,Omega-LabcOmega==0];
%             F = [F,W+Omega-Wcc-Wdc==0];
% 
%             dO = dAO*dBO*dCO;
%             J=real(trace(Omega/dO));
%         elseif isequal(robustness,'RR')
%             J = sdpvar(1,1);
%             F = [F,J>=0,W+J*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI)-Wdc-Wcc==0];
%         end
% 
%         ops = sdpsettings('solver','mosek','verbose',1);
%         ops.mosek.MSK_IPAR_NUM_THREADS=6;
%         SOLUTION=optimize(F,J,ops)
% 
%         varargout{1} = value(J);
%         varargout{2} = value(Wdc);
%         varargout{3} = value(Wcc);
%         if isequal(robustness,'GR')
%             varargout{4} = value(Omega);
%         end
%     end
end