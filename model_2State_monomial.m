% Function returns model of 2 state example in [1] in form of the literature. 


% [Reference]
%   [1] F. Amato, C. Cosentino, and A. Merola, “On the region of asymptotic
% 	stability of nonlinear quadratic systems,” in 2006 14th Mediterranean
% 	Conference on Control and Automation, pp. 1–5, IEEE, 2006.

function model = model_2State_monomial()
	% system dimension
	nx = 2;
	nz = 1;

	model.nx = nx;
	model.nz = nz;

	%% Dynamics
	model.A = [-50, -16; 13 -9];

    % quadratic part of the dynamics
    model.B = [13.8; 5.5];
    model.func_QuadDyn = @(i) QuadDyn(model, i);
    model.Qis = zeros(nx,nx,nz);
    for i = 1:nz
    	model.Qis(:,:,i) = model.func_QuadDyn(i);
    end

    % ODE of the model
    model.func_ode = @(t,x) ode_2State(model,t,x);
    model.func_ode_backward = @(t,x) ode_2State(model,-t,x);

    %% functions to get Quadratic constraints 
    % Meq:	equality QC
    model.Meq = [];     % doesn't have any equality constraints
    
    % Mi: 	inequality QC
    model.Mi.CS_z = @(E) func_getQCs(model,'Mi_CS_z',E);
    model.Mi.Valley_z = @(E) func_getQCs(model,'Mi_Valley_z',E);
end

%% return i-th row of quadratic term of dynamics
%   xdot(i) = A(i,:)*x + x'*Qi*x
function Qi = QuadDyn(model,i)
    nx = model.nx;
    
    % initialize Qi
    Qi = zeros(nx);
    
    switch i
        case 1
            Qi = [0 1; 1 0]/2;
    
        otherwise
            error('Out of system dimension');    
    end
end

%% return ode of the system
function xdot = ode_2State(model,t,x)
    nz = model.nz;
    
    % initialize xdot
    Qterms = zeros(nz,1);
    
    % store the quadratic terms
    for i = 1:nz
        Qterms(i) = x'*model.Qis(:,:,i)*x;
    end
    
    xdot = sign(t)* (model.A*x + model.B*Qterms);    % allow simulation backward in time
end