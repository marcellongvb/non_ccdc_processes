% SDP to find optimal witness for the A->B/B->A scenario
% for process matrices
% 
% Inputs:
% 
% d = vector of dimensions [dAI dAO dBI dBO dCI dCO], for bipartite
% scenario, use dCI = dCO = 1
% 
% robustness = 'GR' for generalized robustness or 'RR' for random
% robustness
% 
% sdp_prob = 'primal' for finding optimal noise Omega such that W + Omega =
% Separable. 'dual' for finding optimal witness operator S that minimizes
% the inner product tr(S*W)
% 
% W: Input process matrix

function [varargout] = order_sdp(d,robustness,sdp_prob,W)

% Input dimensions
dAI = d(1);
dAO = d(2);
dBI = d(3);
dBO = d(4);
dCI = d(5);
dCO = d(6);

if nargin < 4    
    if dCI > 1
    %   Tripartite Scenario, If W is not given, the SWITCH is used as
    %   standard example
    disp('Tripartite Scenario, process SWITCH is used')
    zero = [1; 0];
    one = [0; 1];

    dCIt = 2;
    dCIc = 2;

    psi = 1/(sqrt(2)) * (zero + one);
    phiplus = kron(zero,zero) + kron(one, one);

    term_1 = kron2(psi,phiplus,phiplus,zero);

    term_2 = kron2(psi,phiplus,phiplus,one);
    term_2 = PermuteSystems(term_2,[3 4 1 2 5 6],[dBI dBO dAI dAO dCIt dCIc]);
    %    W Switch:
    W = 1/sqrt(2) * (term_1 + term_2);
    W = W * W';
    clear term_1 term_2 psi phiplus zero one
    elseif dCI == 1
    % Bipartite Scenario, If W is not given, Wocb is used as standard
    % example
        disp('Bipartite Scenario, Wocb is used')
        Z = [1, 0; 0, -1];
        X = [0, 1; 1, 0];

        W = 1/4 * (eye(dAI*dAO*dBI*dBO)+ 1/sqrt(2) * (kron2(eye(2),Z,Z,eye(2))+ kron2(Z,eye(2),X,Z)));
    end
end

if isequal(sdp_prob,'primal')
    disp('sdp_prob set as PRIMAL')
    disp('Finding Witness of Causal Nonseparability')
    yalmip('clear');
    S  = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
    S1 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
    S2 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
    S3 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
    S4 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');

    % Labc(W) = CO - CICO + BOCICO - BIBOCICO + AOBIBOCICO

    LabcS3 = traceandrep(S3,6,d) - traceandrep(S3,[5 6],d) + traceandrep(S3,[4 5 6],d) ...
    - traceandrep(S3,[3 4 5 6],d) + traceandrep(S3,[2 3 4 5 6],d);

    % Lbac(W) = CO - CICO + AOCICO - AIAOCICO + AIAOBOCICO

    LbacS4 = traceandrep(S4,6,d) - traceandrep(S4,[5 6],d) + traceandrep(S4,[2 5 6],d) ...
    - traceandrep(S4,[1 2 5 6],d) + traceandrep(S4,[1 2 4 5 6],d);

    %  Constraints:
    F = [S1>=0,S2>=0,S-S3-S1==0,S-S2-S4==0,LabcS3==0,LbacS4==0];

    if isequal(robustness,'GR') % Generalized Robustness
        disp('Witness: Generalized Robustness')

        S5 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
        S6 = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');

        % AO - AI AO + BO - AO BO + AI AO BO - BI BO + AO BI BO - 
        %  AI AO BI BO + CO - AO CO + AI AO CO - BO CO + AO BO CO - 
        %  AI AO BO CO + BI BO CO - AO BI BO CO + AI AO BI BO CO - CI CO + 
        %  AO CI CO - AI AO CI CO + BO CI CO - AO BO CI CO + AI AO BO CI CO - 
        %  BI BO CI CO + AO BI BO CI CO
        
        if dCI > 1
        LvS5 = validprocessproj(S5,d);
        elseif dCI == 1
           LvS5 = traceandrep(S5,2,d) + traceandrep(S5,4,d) - traceandrep(S5,[2 4],d) - traceandrep(S5,[3 4],d) + traceandrep(S5,[2 3 4],d) ...
        - traceandrep(S5,[1 2],d) + traceandrep(S5,[1 2 4],d);
        end
        
        F = [F,S6>=0,LvS5==0,(eye(dAI*dAO*dBI*dBO*dCI*dCO)/(dAO*dBO*dCO))-S-S5-S6==0];

    elseif isequal(robustness,'RR') % Random Robustness
        disp('Witness: Random Robustness')
        F = [F,trace(S)<=dAI*dBI*dCI]; 
    end

    J=real(trace(S*W));

    ops = sdpsettings('solver','mosek','verbose',1);
    ops.mosek.MSK_IPAR_NUM_THREADS=6;
    SOLUTION=optimize(F,J,ops)

    varargout{1} = value(J);
    varargout{2} = value(S);

elseif isequal(sdp_prob,'dual')
    disp('sdp_prob set as PRIMAL')
    disp('Finding processes to add to the given process and turn it into a separable process')
    
    yalmip('clear');    
    Wab   = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');
    Wba   = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dBI*dCI*dAO*dBO*dCO, 'hermitian', 'complex');

    %  Projector onto linear subspace of Process Matrices compatible with
    %  A->B->C:
    Lab = traceandrep(Wab,6,d) - traceandrep(Wab,[5 6],d) + traceandrep(Wab,[4 5 6],d)- traceandrep(Wab,[3 4 5 6],d) + traceandrep(Wab,[2 3 4 5 6],d);

    %  Projector onto linear subspace of Process Matrices compatible with
    %  B->A->C:

    Lba = traceandrep(Wba,6,d) - traceandrep(Wba,[5 6],d) + traceandrep(Wba,[2 5 6],d) - traceandrep(Wba,[1 2 5 6],d) + traceandrep(Wba,[1 2 4 5 6],d);
    
    F = [Wab>=0,Wba>=0,Lab-Wab==0,Lba-Wba==0];

    if isequal(robustness,'GR')
        disp('Generalized Robustness')
        Omega = sdpvar(dAI*dAO*dBI*dBO*dCI*dCO,dAI*dAO*dBI*dBO*dCI*dCO, 'hermitian', 'complex');

        if dCI > 1
            Lv = validprocessproj(Omega,d);
        elseif dCI == 1
            Lv = traceandrep(Omega,2,d) + traceandrep(Omega,4,d) - traceandrep(Omega,[2 4],d) - traceandrep(Omega,[3 4],d) + traceandrep(Omega,[2 3 4],d) ...
            - traceandrep(Omega,[1 2],d) + traceandrep(Omega,[1 2 4],d);
        end
            
        F = [F,W+Omega-Wab-Wba==0,Omega>=0,Lv-Omega==0];

        % Objective function
        J = real(trace(Omega/(dAO*dBO)));
   
    elseif isequal(robustness,'RR')
        disp('Random Robustness')
        J = sdpvar(1,1);
        
        F = [F,J>=0,W+J*eye(dAI*dAO*dBI*dBO*dCI*dCO)/(dAI*dBI*dCI)-Wab-Wba==0];
 
    end
        ops = sdpsettings('solver','mosek','verbose',1);
        % ops.mosek.MSK_IPAR_NUM_THREADS=6;
        SOLUTION = optimize(F,J,ops)

        varargout{1} = value(J);        
        varargout{2} = value(Wab);
        varargout{3} = value(Wba);
        
        if isequal(robustness,'GR')
            varargout{4} = value(Omega);
        end
end
end