function [R] = ccdc_separability_outer(W, dims, sdp_prob, robustness, k)
% CCDC_SEPARABILITY MEASURED USING OUTER APPROXIMATION OF CCDC SET
% This function is used for calculating the CCDC robustness of a given
% process W with dimensions (dAI, dAO, dBI).
% Arguments:
% W -- Input process 
% dims -- Vector of dimensions [dAI dAO dBI]
% sdp_prob -- options: 'primal','dual'. 'primal' for calculating robustness 
% via decomposition problem, 'dual'for calculating robustness via witness 
% of CCDC separability.
% robustness -- options: 'GR', 'WNR'. 'GR' for generalised robustness,
% 'WNR' for white noise robustness.
% k -- level of symmetric extensions

% Dimensions of subspaces
dAI = dims(1);
dAO = dims(2);
dBI = dims(3);

yalmip('clear');
if isequal(sdp_prob,'primal')   
%     Common-cause term:
        rho = sdpvar(dAI*dBI, dAI*dBI, 'hermitian', 'complex');
        Wcc = PermuteSystems(Tensor(rho,eye(dAO)),[1 3 2], [dAI dBI dAO]);
        F = [rho>=0];
%         Direct-cause term:
        Wdc = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');
        Wdc_k = sdpvar((dAI^k)*dAO*dBI, (dAI^k)*dAO *dBI, 'hermitian', 'complex');
        Projdc = Wdc - traceandrep(Wdc, 3, [dAI dAO dBI]) + traceandrep(Wdc, [2 3], [dAI dAO dBI]);
        F = [F,Wdc>=0, Wdc-Projdc==0, PartialTrace(Wdc_k, 2, [dAI dAI^(k-1) dAO dBI])-Wdc==0, PartialTranspose(Wdc_k, 1, [dAI dAI^(k-1) dAO dBI])>=0];
        F = [F, Tensor(full(SymmetricProjection(dAI,k)),eye(dAO*dBI))*Wdc_k == Wdc_k]; %Bose symmetric extension
    if isequal(robustness, 'GR')
%         Generalised Noise
        Omega = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');
        F = [F, Omega>=0,(1-trace(Omega)/dAO)*W+Omega-Wcc-Wdc==0];
%         Objective
        R = real(trace(Omega)/dAO);        
        ops = sdpsettings('solver', 'mosek', 'verbose', 1);
        ops.mosek.MSK_IPAR_NUM_THREADS=6;
        SOLUTION=optimize(F,R,ops);        
        R = value(R);              
    elseif isequal(robustness, 'WNR')
        R = sdpvar(1,1,'hermitian','complex');
        F = [F,(1-R)*W+R*eye(dAI*dAO*dBI)/(dAI*dBI)-Wcc-Wdc==0];        
        ops = sdpsettings('solver', 'mosek', 'verbose', 1);
        ops.mosek.MSK_IPAR_NUM_THREADS=6;
        SOLUTION=optimize(F,R,ops);        
        R = value(R);           
    end        
elseif isequal(sdp_prob, 'dual')   
%     For the dual problems, we only use the most simple case k=1 for
%     implementing the see-saw
    S = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');
    Sdc = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');
    Sdcab = sdpvar(dAI*dAO*dBI, dAI*dAO*dBI, 'hermitian', 'complex');    
    Projab = Sdcab + traceandrep(Sdcab, [2 3], dims) - traceandrep(Sdcab, 3, dims);
    if isequal(robustness, 'GR')          
        F = [Sdc>=0,S-PartialTranspose(Sdc, 1, dims)+Sdcab-Projab>=0, PartialTrace(S,2,[dAI dAO dBI])>=0];   
        F = [F, eye(dAI*dAO*dBI)*(1+trace(S*W))-dAO*S>=0];        
    elseif isequal(robustness, 'WNR')
        F = [trace(S)-(dAI*dBI)*trace(S*W)<=(dAI*dBI),PartialTrace(S,2,[dAI dAO dBI])>=0,S-PartialTranspose(Sdc,1,[dAI dAO dBI])-Sdcab+Projab>=0,Sdc>=0];        
    end
    R = trace(S*W);        
    ops = sdpsettings('solver', 'mosek', 'verbose', 1);
    ops.mosek.MSK_IPAR_NUM_THREADS=6;
    SOLUTION=optimize(F,R,ops);        
    R = -value(R);  
    
end
end



