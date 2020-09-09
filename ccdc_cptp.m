function [varargout] = ccdc_cptp(d,p,robustness,W)

 % Input dimensions
    dAI = d(1);
    dAO = d(2);
    dBI = d(3);
    dBO = d(4);
    dCI = d(5);
    dCO = d(6);

    J = p;
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

        F = [F,Omega>=0,Omega-LabcOmega==0,trace(Omega)==dAO*dBO*dCO];
        F = [F,(1-p)*W+p*Omega==Wcc+Wdc];

        dO = dAO*dBO*dCO;
%         J=real(trace(Omega/dO));
    elseif isequal(robustness,'RR')
%         J = sdpvar(1,1);
%         F = [F,J>=0,(1-p)*W+p*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI)==Wdc+Wcc];
        F = [F,(1-p)*W+p*eye(dAI*dAO*dBI*dBO*dCI)/(dAI*dBI*dCI)==Wdc+Wcc];
    end

    ops = sdpsettings('solver','mosek','verbose',1);
    ops.mosek.MSK_IPAR_NUM_THREADS=6;
    SOLUTION=optimize(F,J,ops)

    varargout{1} = value(J);
    varargout{2} = value(Wdc);
    varargout{3} = value(Wcc);
    if isequal(robustness,'GR')
        varargout{4} = value(Omega);
    end


end