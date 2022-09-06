


function [r, info] = SDPOpt_ROA_Etile(model, xi, Meq, Mi, opt)
    %% check options
    if nargin < 5
        opt = func_getOptions_SDP_ROA();
    end

    %% setup
	% system dynamics
	A = model.A;
	B = model.B;
	nx = model.nx;
	nz = model.nz;

	% utility matrices
	Inx = eye(nx);
	Znx = zeros(nx);
	Znx1 = zeros(nx,1);

	Inz = eye(nz);
    Znz = zeros(nz);
    Znz1 = zeros(nz,1);

    % find number of local quadratic constraints
    if isempty(Meq)
        NMeq = 0;
    else
        NMeq = size(Meq,3);
    end
    if isempty(Mi)
        NMi = 0;
    else
        NMi = size(Mi.tilde,3);
        
        Mi_tilde = Mi.tilde;
        Mi_hat = Mi.hat;
        
        Mi_L = Mi.L;
    end
    
    %% solve SDP with Etilde fixed
    cvx_begin sdp quiet
    cvx_solver mosek
        variable P(nx,nx)       semidefinite;
        variable EtildeInv(nx,nx)  semidefinite;
        variable nu(NMeq,1);
        variable lam            nonnegative;
        
        minimize(lam);
        
        subject to 
            % set containment;
            P <= lam*Inx;
            [P      Inx;
             Inx    EtildeInv] >= 0;
            
            % Lyapunov stability
            Spro = [A'*P+P*A+opt.eps*Inx    P*B;
                    B'*P                    Znz];
            for i = 1:NMeq
                Spro = Spro + nu(i)*Meq(:,:,i);
            end
            for i = 1:NMi
                L = Mi_L(:,:,i);
                Mi_tilde = blkdiag(L*EtildeInv*L', Znz);
                
                Spro = Spro + xi(i)*(Mi_tilde + Mi_hat(:,:,i));
            end
            Spro <= 0;
    cvx_end
    
    % set output
    info = get_infoStru();
    info.xi = xi;
    
    % check feasibility
    if strcmp(cvx_status, 'Solved')
        % solve for r
        r = 1/sqrt(lam);
        if opt.verbose
            fprintf('[E-step]  r = %2.6f \n', r);
        end
        
        % save info
        info.feasibility = true;
        info.r = r;
        info.Etilde = inv(EtildeInv);
        info.P = P;
        info.nu = nu;
        info.LMI = Spro;     
    else
        r = -Inf;
        if opt.verbose
            fprintf('[E-step]  Infeasible!');
        end
        
        % save info
        info.feasibility = false;
        info.r = r;
    end
end