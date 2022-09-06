






function model = model_WKH(params)
    if nargin < 1
       params = struct('ParamSet','B&T','lambda',1,'mu',1,'nu',1,'rho',1,'delta',1,'gamma',1);
       params.Re = 1;       
    end

    nx = 4;
    nz = 5;
    
    model.nx = nx;
    model.nz = nz;
    
    %% Dynamics
    model.params = params;
    
    % linear part of WKH dynamics
    A = -diag([params.lambda,params.mu,params.nu,params.rho])/params.Re;
    A(1,2) = 1;
    model.A = A;
    
    % Quadratic part of WKH dynamics
    % z = [x3^2, x1x2, x1x3, x2x3, x2x4] 
    B = zeros(nx,nz);
    B(1,1) = -params.gamma;
    B(1,5) = 1;
    B(2,1) = params.delta;
    B(3,3) = params.gamma;
    B(3,4) = -params.delta;
    B(4,2) = -1;
    model.B = B;
    
    % quadratic part of WKH dynamics
    model.func_QuadDyn = @(i) QuadDyn(model,i);
    model.Qis = zeros(nx,nx,nz);
    for i = 1:nz
        model.Qis(:,:,i) = model.func_QuadDyn(i);
    end
    
    model.Qhats = zeros(nx,nx,nx);
    for i = 1:nx
        model.Qhats(:,:,i) = get_Qhat(model,i);
    end
    
    % ODE of WKH dynamics
    model.dyn = @(t,x) WKH_dyn(model,t,x);
    
    %% qudratic constraints
    % [Meq] 
    model.Meq.lossless = @(E) func_getQCs(model,'Meq_lossless',E);
    
    % [Mi]
    model.Mi.CS_z = @(E) func_getQCs(model, 'Mi_CS_z',E);
    
%     CS_phi1 = @(E) getQC_CS_phi1(model,E);
%     CS_phi3 = @(E) getQC_CS_phi3(model,E);
%     tilde = @(E) cat(3, CS_phi1(E).tilde, CS_phi3(E).tilde);
%     hat = @(E) cat(3, CS_phi1(E).hat, CS_phi3(E).hat);
%     model.Mi.CS_phi = @(E) struct('tilde',tilde(E), 'hat',hat(E));
    model.Mi.CS_phi = @(E) func_getQCs(model, 'Mi_CS_phi',E);
    
    model.Mi.Valley2_z = @(E) func_getQCs(model, 'Mi_Valley_z',E);
    
%     model.Mi.Valley2_phi = @(E) getQC_VL2_phi(model,E);
    model.Mi.Valley2_phi = @(E) func_getQCs(model, 'Mi_Valley_phi',E);
    
    model.Mi.Valley3_phi = @(E) getQC_VL3(model,E);
    
    model.Mi.CrossP_z = @(E) getQC_CrossP_z(model,E);
%     model.Mi.CrossP_neg = @(E) getQC_CrossP_z_negative(model,E);
end

%% return i-th row of quadratic term of dynamics
%   xdot(i) = A(i,:)*x + x'*Qi*x
function Qi = QuadDyn(model,i)
    nx = model.nx;
    
    % initialize Qi
    Qi = zeros(nx);
    
    switch i
        case 1      % z1 = x3^2
            Qi(3,3) = 1;
            
        case 2      % z2 = x1*x2
            Qi(1,2) = 1/2;
            Qi(2,1) = 1/2;
            
        case 3      % z3 = x1*x3
            Qi(1,3) = 1/2;
            Qi(3,1) = 1/2;
            
        case 4      % z4 = x2*x3
            Qi(2,3) = 1/2;
            Qi(3,2) = 1/2;
            
        case 5      % z5 = x2*x4
            Qi(2,4) = 1/2;
            Qi(4,2) = 1/2;
            
        otherwise
            error('quadratic dynamic not implemented'); 
    end
end

function Qhat = get_Qhat(model,i)
    nx = model.nx;
    B = model.B;
    
    % initialize Qhat
    Qhat = zeros(nx);
    
    for j = 1:model.nz
        Qhat = Qhat + B(i,j)*model.Qis(:,:,j);
    end
end

function xdot = WKH_dyn(model,t,x)
    % notation follow eq(2) in [1]
    
    % system dimension
    nz = model.nz;

    % nonlinearity (quadratic part) of the dynamics
    z = zeros(nz,1);
    for i = 1:nz
        z(i) = x'*model.Qis(:,:,i)*x;
    end
    
    xdot = model.A*x + model.B*z;    
end

% ### Hand-craft
function Mi = getQC_CS_phi3(model,E)
    % z3 = x1*x3;
    % z4 = x2*x3;

    % nonlinearity in 3rd row: 
    %   y3 = gamma*z3 - delta*z4
    
    % Cauchy-Schwarz for y3:
    %   y3^2 = (gamma*z3 - delta*z4)^2                          ---*
    %        = x'* (gamma*Q3 - delta*Q4) *x
    %        <= ||x||_2^2 * ||(gamma*Q3 - delta*Q4) *x||_2^2
    %        <= R^2 *x'* (gamma*Q3 - delta*Q4)^2 *x             ---*
    
    % Form the QC in terms of [x; z]
    
    
    % dimension of the system
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % property of system
%     gamma = model.params.gamma;
%     delta = model.params.delta;
%     Q3 = model.func_QuadDyn(3);
%     Q4 = model.func_QuadDyn(4);
    
    B = model.B;
    
    Qis = model.Qis;

    % initialization 
    Mi.tilde = zeros(nxz);
    Mi.hat = zeros(nxz);
    
    % forming the QC    
%     Mi.tilde(1:nx,1:nx) = (gamma*Q3 - delta*Q4)^2;
    % ########## this is probably sigma( B(j,:).*Qis(:,:,i))
    % yj to be bound
    % zj are qis
    Qy3 = Znx;
    for i = 1:nz
        Qy3 = Qy3 + B(3,i)*Qis(:,:,i);
    end
    Mi.tilde(1:nx,1:nx) = Qy3*inv(E)*Qy3;
    
%     Mi.hat(nx+3,nx+3) = -gamma^2;
%     Mi.hat(nx+3,nx+4) = gamma*delta;
%     Mi.hat(nx+4,nx+3) = gamma*delta;
%     Mi.hat(nx+4,nx+4) = -delta^2;
    % ########## this is prabably -B'*B for the particular row!
    B3 = B(3,:);
    Mi.hat(nx+1:end,nx+1:end) = -B3'*B3;
end

% ### Hand-craft
function Mi = getQC_CS_phi1(model,E)
    % dimension of the system
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % property of system 
    B = model.B;
    Qis = model.Qis;

    % initialization 
    Mi.tilde = zeros(nxz);
    Mi.hat = zeros(nxz);
    
    % forming the QC    
%     Mi.tilde(1:nx,1:nx) = (gamma*Q3 - delta*Q4)^2;
    % ########## this is probably sigma( B(j,:).*Qis(:,:,i))
    % yj to be bound
    % zj are qis
    Qy1 = Znx;
    for i = 1:nz
        Qy1 = Qy1 + B(1,i)*Qis(:,:,i);
    end
    Mi.tilde(1:nx,1:nx) = Qy1*inv(E)*Qy1;
    
%     Mi.hat(nx+3,nx+3) = -gamma^2;
%     Mi.hat(nx+3,nx+4) = gamma*delta;
%     Mi.hat(nx+4,nx+3) = gamma*delta;
%     Mi.hat(nx+4,nx+4) = -delta^2;
    % ########## this is prabably -B'*B for the particular row!
    B1 = B(1,:);
    Mi.hat(nx+1:end,nx+1:end) = -B1'*B1;
end

% function lmax = func_lambdaMax(Q)
% % compute the largest eigenvalues of matrix Q
%     lmax = max(eig(Q));
% 
%     if lmax <= 0
%         warning('No positive eigenvalue!');
%     end
% end


% ### Hand-craft
function Mi = getQC_VL3(model,E)
% rank 3 valley QC on phi:      phi_1 = x2*x4 - x_3^2
%   perform QC on -phi_1!

    % precompute matrices
    invE = inv(E);
    invEsq = sqrtm(invE);

    % dimension of the system
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % initialization of output
    Mi.tilde = zeros(nxz,nxz,4);
    Mi.hat = zeros(nxz,nxz,4);
    W = zeros(nx,nx,4);
       
    % compute Qp for phi = x'Qp x = bZ
    Qp = Znx;
    for i = 1:nz
        Qp = Qp + model.B(1,i)*model.Qis(:,:,i);
    end
    
%     Q = [0 0 0 0; 0 0 0 0.5; 0 0 -1 0; 0 0.5 0 0];

    % eigenvalues of -Qp
    eig1 = 0.5;
    eig2 = 1;
    eig3 = -0.5;
    
    v1 = [0 -1/sqrt(2) 0 1/sqrt(2)]';
    v2 = [0 0 1 0]';
    v3 = [0 -1/sqrt(2) 0 -1/sqrt(2)]';
    
    % compute axillary vectors
    nu1 = sqrt(eig1)*v1;
    nu2 = sqrt(eig2)*v2;
    nu3 = sqrt(abs(eig3))*v3;
    
    % first set of QC
    c1 = nu1 + nu3;
    c2 = nu1 - nu3;
    S = nu2*nu2';
    
    beta1 = c2'*invE*c2;
    beta2 = c1'*invE*c1;
    gamma = func_lambdaMax(invEsq*(2*Qp - S)*invEsq);
    
    W(:,:,1) = beta1*(c1*c1') + gamma*S;
    W(:,:,2) = beta2*(c2*c2') + gamma*S;
    
    % second set of QC
    c1 = nu2 + nu3;
    c2 = nu2 - nu3;
    S = nu1*nu1';
    
    beta1 = c2'*invE*c2;
    beta2 = c1'*invE*c1;
    gamma = func_lambdaMax(invEsq*(2*Qp - S)*invEsq);
    
    W(:,:,3) = beta1*(c1*c1') + gamma*S;
    W(:,:,4) = beta2*(c2*c2') + gamma*S;
    
    % setting output
    B1 = model.B(1,:);
    
    for i = 1:4
        Mi.tilde(:,:,i) = blkdiag(W(:,:,i), Znz);
        Mi.hat(:,:,i) = blkdiag(Znx, -B1'*B1);              %%% check this!!!
    end
end

% ### Hand-craft
function Mi = getQC_VL2_phi(model,E)
% rank 2 valley QC on phi:      phi_3 = x1*x3 - x2*x3
    
    % system parameter
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % initialization of output
    Mi.tilde = zeros(nxz,nxz,2);
    Mi.hat = zeros(nxz,nxz,2);
    
    % forming Qhat3 and b3
    b3 = model.B(3,:);
    
    Qhat3 = Znx;
    for i = 1:nz
        Qhat3 = Qhat3 + b3(i)*model.Qis(:,:,i);
    end    

    
    % get W from func_getW_Valley
    W = func_getW_Valley2(model, Qhat3, E);
    
    % setting output
    for i = 1:2
        Mi.tilde(:,:,i) = blkdiag(W(:,:,i), Znz);
        Mi.hat(:,:,i) = blkdiag(Znx, -b3'*b3);
    end
end

% ### Hand-craft
function Mi = getQC_CrossP_z(model,E)
% cross product on Z = [x3^2, x1x2, x1x3, x2x3, x2x4] 
% - z1*zi                       #4
% - z2*z3 = x1^2 * x2 * x3
% - z2*z5 = x1 * x2^2 * x4
% - z3*z4 = x1 * x2 * x3^2
% - z4*z5 = x2^2 * x3 * x4

    % system dimension
    nx = model.nx;
    nz = model.nz;
    
    nxz = nx + nz;
    
    % initialize output
    Mi.tilde = zeros(nxz,nxz,8*2);
    Mi.hat = zeros(nxz,nxz,8*2);

    % setting vectors for func_getCrossProductQC
    p_list = [1, 1, 1, 1, 2, 2, 3, 4];
    q_list = [2, 3, 4, 5, 3, 5, 4, 5];
    i_list = [3, 3, 3, 3, 1, 2, 3, 2];
    j_list = [1, 1, 2, 2, 2, 1, 1, 3];
    k_list = [2, 3, 3, 4, 3, 4, 2, 4];
    
    for idx = 1:8
        p = p_list(idx);
        q = q_list(idx);
        i = i_list(idx);
        j = j_list(idx);
        k = k_list(idx);
        
        M = func_getCrossProductQC(model,E, p,q,i,j,k);
        
        idM = 4*(idx-1) + [1,2,3,4];
        Mi.tilde(:,:,idM) = M.tilde;
        Mi.hat(:,:,idM) = M.hat;
    end
end

function Mi = getQC_CrossP_z_negative(model,E)
% cross product on Z = [x3^2, x1x2, x1x3, x2x3, x2x4] 
% - z1*zi                       #4
% - z2*z3 = x1^2 * x2 * x3
% - z2*z5 = x1 * x2^2 * x4
% - z3*z4 = x1 * x2 * x3^2
% - z4*z5 = x2^2 * x3 * x4

    % system dimension
    nx = model.nx;
    nz = model.nz;
    
    nxz = nx + nz;
    
    % initialize output
    Mi.tilde = zeros(nxz,nxz,8*2);
    Mi.hat = zeros(nxz,nxz,8*2);

    % setting vectors for func_getCrossProductQC
    p_list = [1, 1, 1, 1, 2, 2, 3, 4];
    q_list = [2, 3, 4, 5, 3, 5, 4, 5];
    i_list = [3, 3, 3, 3, 1, 2, 3, 2];
    j_list = [1, 1, 2, 2, 2, 1, 1, 3];
    k_list = [2, 3, 3, 4, 3, 4, 2, 4];
    
    for idx = 1:8
        p = p_list(idx);
        q = q_list(idx);
        i = i_list(idx);
        j = j_list(idx);
        k = k_list(idx);
        
        M = func_getCrossProductQC_negative(model,E, p,q,i,j,k);
        
        idx_M = 2*(idx-1) + [1,2];
        Mi.tilde(:,:,idx_M) = M.tilde;
        Mi.hat(:,:,idx_M) = M.hat;
    end
end