function [Lp_hv1, Lq_hv1] = WaterFill_KF(L, P, Q, Param)
    tp = Param.tp;
    tq = Param.tq;
    
    % Duplication matrices
    DmP = duplication_matrix(P); % NC
    DmQ = duplication_matrix(Q);
    
    % Constraints
    Cm1 = [vec(eye(P))'*DmP; kron(ones(P,1)',eye(P))*DmP];  % NC
    Cm2 = [vec(eye(Q))'*DmQ; kron(ones(Q,1)',eye(Q))*DmQ];  % NC
    C = [Cm1, zeros(P+1, 0.5*Q*(Q+1)); % NC
        zeros(Q+1,0.5*P*(P+1)), Cm2];
    
    dp = [tp ;zeros(P,1)];
    dq = [tq ;zeros(Q,1)];
    d = [tp ;zeros(P,1); tq ;zeros(Q,1)];
    Hm = 2*blkdiag( Q*DmP'*DmP,  P*DmQ'*DmQ);
    
    Ltilde = TildeTransform(L,Q,Q,P,P);
    qv = -2*[vec(eye(Q))'*Ltilde*DmP, vec(eye(P))'*Ltilde'*DmQ]';
    
    v =0.5*P*(P+1) + 0.5*Q*(Q+1); % Number of variable we will solve for QP
    % cvx_begin quiet
    % variable z(v,1)
    % minimize ( (1/2)*quad_form(z,Hm) + qv'*z )
    % subject to
    % C*z == d
    % z >= 0
    % cvx_end
    %
    % lp_c_hat = z(1:0.5*P*(P+1));
    % lq_c_hat = z((0.5*P*(P+1))+1:end);
    %
    % Lp_hv1 = reshape(DmP*lp_c_hat,P,P);
    % Lq_hv1 = reshape(DmQ*lq_c_hat,Q,Q);
    
    
    %% Waterfilling based solver
    
    [z1 err] = PGL(Hm, qv, C, d, 1e-6, 0.1);
    lp_c_hat1 = z1(1:0.5*P*(P+1));
    lq_c_hat1 = z1((0.5*P*(P+1))+1:end);
    
    Lp_hv1 = full(reshape(DmP*lp_c_hat1,P,P));  Lp_hv1(abs(Lp_hv1)<10^-4) = 0;
    Lq_hv1 = full(reshape(DmQ*lq_c_hat1,Q,Q));  Lq_hv1(abs(Lq_hv1)<10^-4) = 0;
    
    end
