

function model = model_CDC3States()

    % system dimension
    nx = 3;
    nz = nx*(nx+1)/2;
    
    model.nx = nx;
    model.nz = nz;
    
    %% Dynamics
    % xdot = Ax + Bw,   w is quadratic monomials
    model.A = [-1 -1 -1; -1 -6 -1; -1 -1 -13];
    
    B = zeros(nx,nz);
    B(1,2) = 1;
    B(1,6) = 1;
    B(2,3) = -6;
    B(2,5) = -4;
    B(3,5) = 2;
    B(3,6) = -1;
    model.B = B;
    
    % quadratic monomials
    % Z(x) = [x1^2, x1x2, x1x3, x2^2, x2x3, x3^2];
    model.func_QuadDyn = @(i) QuadDyn(model,i);
    model.Qis = zeros(nx,nx,nz);
    for i = 1:nz
        model.Qis(:,:,i) = model.func_QuadDyn(i);
    end
    
    model.Qhats = zeros(nx,nx,nx);
    for i = 1:nx
        model.Qhats(:,:,i) = get_Qhat(model,i);
    end
    
    % ODE of the dynamics
    model.dyn = @(t,x) func_ODE(model,t,x);
    
    %% Quadratic constraints
    % [Meq]
    model.Meq = [];
    
    % [Mi]
    % CSQC
    model.Mi.CS_z = @(E) func_getQCs(model, 'Mi_CS_z', E);
    model.Mi.CS_phi = @(E) func_getQCs(model, 'Mi_CS_phi', E);
    
    % Valley QC rank2
    model.Mi.Valley2_z = @(E) func_getQCs(model, 'Mi_Valley_z', E);
    model.Mi.Valley2_phi = @(E) func_getQCs(model, 'Mi_Valley_phi', E);
    
    % Valley QC rank3
    model.Mi.Valley3_phi = @(E) func_getQCs(model, 'Mi_Valley3_phi', E);
    
    % Cross-product QC
    model.Mi.CrossP_z = @(E) func_getQCs(model, 'Mi_CrossProduct_z', E);    
%     model.Mi.CrossP_neg = @(E) func_getQCs(model, 'Mi_CrossProduct_neg', E);
%     model.Mi.CrossP_xixj = @(E) func_getQCs(model, 'Mi_CrossProduct_xixj', E);
end

%% Dynamics
function Qi = QuadDyn(model,i)
% Z(x) = [x1^2, x1x2, x1x3, x2^2, x2x3, x3^2];
    nx = model.nx;
    
    % initialize Qi
    Qi = zeros(nx);
    
    % base on each monomial
    switch i
        case 1  % x1^2
            Qi(1,1) = 1;
            
        case 2  % x1x2
            Qi(1,2) = 0.5;
            Qi(2,1) = 0.5;
            
        case 3  % x1x3
            Qi(1,3) = 0.5;
            Qi(3,1) = 0.5;
            
        case 4  % x2^2
            Qi(2,2) = 1;
            
        case 5  % x2x3
            Qi(2,3) = 0.5;
            Qi(3,2) = 0.5;
            
        case 6  % x3^2
            Qi(3,3) = 1;
            
        otherwise
            error('quadratic dynamic not implemented');
    end
end

function Qhat = get_Qhat(model,i)
% get quadratic form of nonlinearity in i-th state derivatives
% xdot = Ax + Bz

    % system properties
    nx = model.nx;
    B = model.B;
    
    % initialize Qhati
    Qhat = zeros(nx);
    
    for j = 1:model.nz
        Qhat = Qhat + B(i,j)*model.Qis(:,:,j);
    end
end

function xdot = func_ODE(model,t,x)
    nz = model.nz;
    
    if isnumeric(x)
        w = zeros(nz,1);
    else
        w = sym('w',[nz,1]);
    end
    
    for i = 1:nz
        w(i) = x'*model.Qis(:,:,i)*x;
    end
    
    xdot = model.A*x + model.B*w;
end

%% Quadratic Constraints
% Z(x) = [x1^2, x1x2, x1x3, x2^2, x2x3, x3^2];
%   - w1*w2     1
%   - w1*w3     2
%   - w1*w5     3
%   - w2*w3     4
%   - w2*w4     5
%   - w2*w5     6
%   - w2*w6     7
%   - w3*w4     8
%   - w3*w5     9
%   - w3*w6     10
%   - w4*w5     11
%   - w5*w6     12