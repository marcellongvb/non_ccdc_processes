function [R] = W_FB_robustness_lower(robustness)

% SDP for calculating lower-bound for CCDC robustnesses of the process W_FB
% input d: vector of dimensions
% 
% input robustness: 'GR' for generalized robustness, 'RR' for random
% robustness
% 
% input sdp_prob: 'primal' (for witness) or 'dual' (for decomposition)
% 
% opt input W: Process matrix to be tested. default: Wcoh

    % Input dimensions
    dAI = 2;
    dAO = 2;
    dBI = 2;
    dBO = 2;
    dCI = 8;


    % |W_FB>>:
    % 1/sqrt(2) ( |0>CI(1) |psi>AIBI |identity>AOCI(2) |identity>BOCI(3) +
    % |1>CI(1)|psi>AICI(2) |identity>AOBI |identity>BOCI(3)) 

    zero = [1; 0];
    one = [0; 1];

    dCI1 = 2;
    dCI2 = 2;
    dCI3 = 2;
    psi = 1/(sqrt(2)) * (Tensor(zero,zero) + Tensor(one,one));
    phiplus = Tensor(zero,zero) + Tensor(one, one);

    term_1 = Tensor(zero,psi,phiplus,phiplus);
    term_1 = PermuteSystems(term_1,[2 4 3 6 1 5 7],[dCI1 dAI dBI dAO dCI2 dBO dCI3]);

    term_2 = PermuteSystems(Tensor(one,psi,phiplus,phiplus),[2 4 5 6 1 3 7],[dCI1 dAI dCI2 dAO dBI dBO dCI3]);

    %    W_FB:
    W = 1/sqrt(2) * (term_1 + term_2);
    W = W * W';
    clear term_1 term_2 psi phiplus zero one

    yalmip('clear');

    Wcc = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
    Wdc = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');

    Wdc_transp = PartialTranspose(traceandrep(Wdc,[4 5],d),1,d);
    LabcWdc = Wdc - traceandrep(Wdc,5,d) + traceandrep(Wdc,[4 5],d) - traceandrep(Wdc,[3 4 5],d) +  traceandrep(Wdc,[2 3 4 5],d);
    LabcWcc = Wcc - traceandrep(Wcc,5,d) + traceandrep(Wcc,[4 5],d) - traceandrep(Wcc,[3 4 5],d) +  traceandrep(Wcc,[2 3 4 5],d);


    LccWcc = Wcc - traceandrep(Wcc,5,d) + traceandrep(Wcc,[4 5],d) - traceandrep(Wcc,[3 4 5],d) +  traceandrep(Wcc,[2 3 4 5],d) ...
    - traceandrep(Wcc,[4 5],d) + traceandrep(Wcc,[3 4 5],d) -  traceandrep(Wcc,[2 3 4 5],d) + traceandrep(Wcc,[2 4 5],d);


    F = [Wdc>=0,Wcc>=0,Wdc_transp>=0,Wcc-LabcWcc==0,Wdc-LabcWdc==0,Wcc-LccWcc==0];

    if isequal(robustness,'GR')
        Omega = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');

        LabcOmega = Omega - traceandrep(Omega,5,d) + traceandrep(Omega,[4 5],d) - traceandrep(Omega,[3 4 5],d) +  traceandrep(Omega,[2 3 4 5],d);

        F = [F,Omega>=0,Omega-LabcOmega==0];
        F = [F,(1-trace(Omega)/(dAO*dBO))W+Omega-Wcc-Wdc==0];

        dO = dAO*dBO*dCO;
        R=real(trace(Omega/dO));
    elseif isequal(robustness,'RR')
        R = sdpvar(1,1);
        F = [F,J>=0,(1-R)W+J*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI)-Wdc-Wcc==0];
    end

    ops = sdpsettings('solver','mosek','verbose',1);
    ops.mosek.MSK_IPAR_NUM_THREADS=6;
    SOLUTION=optimize(F,R,ops)
    R = value(R);
    end


