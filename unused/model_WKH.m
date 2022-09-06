






function model = model_WKH(params)
    nx = 4;
    nz = nx;
    
    model.nx = nx;
    model.nz = nz;
    
    %% Dynamics
    model.params = params;
    
    % linear part of WKH dynamics
    A = -diag([params.lambda,params.mu,params.nu,params.rho])/params.Re;
    A(1,2) = 1;
    model.A = A;
    
    B = eye(nx);
    model.B = B;
    
    % quadratic part of WKH dynamics
    model.func_QuadDyn = @(i) QuadDyn(model,i);
    model.Qis = zeros(nx,nx,nz);
    for i = 1:nx
        model.Qis(:,:,i) = model.func_QuadDyn(i);
    end
    
    % ODE of WKH dynamics
    model.dyn = @(t,x) WKH_dyn(model,t,x);
    
    %% qudratic constraints
    % [Meq] 
    model.Meq.lossless = @(E) func_getQCs(model,'Meq_lossless',E);
    
    % [Mi]
    model.Mi.CS_z = @(E) func_getQCs(model, 'Mi_CS_z',E);
end

%% return i-th row of quadratic term of dynamics
%   xdot(i) = A(i,:)*x + x'*Qi*x
function Qi = QuadDyn(model,i)
    nx = model.nx;
    params = model.params;
    
    % initialize Qi
    Qi = zeros(nx);
    
    switch i
        case 1      % z1 = x'*Q1*x (R^1)
            Qi(3,3) = -1/2*params.gamma;
            Qi(4,2) = 1/2;
            Qi = Qi + Qi';
            
        case 2      % z2 = x'*Q2*x (R^1)
            Qi(3,3) = params.delta;
            
        case 3      % z3 = x'*Q3*x (R^1)
            Qi(1,3) = 1/2*params.gamma;
            Qi(2,3) = -1/2*params.delta;
            Qi = Qi + Qi';
            
        case 4      % z4 = x'*Q4*x (R^1)
            Qi(1,2) = -1/2;
            Qi = Qi + Qi';
            
        otherwise
            error('quadratic dynamic not implemented'); 
    end
end

function xdot = WKH_dyn(model,t,x)
    % notation follow eq(2) in [1]
    
    % system dimension
    nz = model.nz;
    
    % store the quadratic terms
    z = zeros(nz,1);
    for i = 1:nz
        z(i) = x'*model.Qis(:,:,i)*x;
    end
    
    xdot = model.A*x + model.B*z;
end