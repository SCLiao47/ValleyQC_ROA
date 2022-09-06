


function [r_max, info_max, cnt_iter, list_r, list_info] = SDPOpt_ROA_ALP_bisection(model, E, Meq, Mi, opt)
	%% check options
	if nargin < 5
        opt = func_getOptions_SDP_ROA();
    end
    if ~isfield(opt, 'tol_alp')
        opt.tol_alp = 0.001;
    end
    
%     opt.verbose = true;
    
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
    end
    
    %% solve SDP with bisection on alpha
    alp_min = max(opt.alp_min, 0.001);
    alp_max = opt.alp_max;
    alp_test = 1;

    % initialize info/storage
    list_alp = [];
    list_r = [];
    list_info = {};
    
    
    % iteration
    cnt_iter = 0;
    del_alp = 1;
    
    while abs(del_alp) > opt.tol_alp
        % iteration count
        cnt_iter = cnt_iter + 1;
        
        % information
        info = get_infoStru();
        info.alpha = alp_test;
        info.E = E;        
        list_alp(cnt_iter) = alp_test;
        
        % [SDP]
        % lam = 1/r^2;
        cvx_begin sdp quiet
        cvx_solver mosek
            variable P(nx,nx)       semidefinite;
            variable nu(NMeq,1);
            variable xi(NMi,1)      nonnegative;
            variable lam            nonnegative;
            
            minimize(lam);
            
            subject to 
                % set containment;
                P <= lam*Inx;
                E/(alp_test)^2 <= P;
                
                % Lyapunov stability condition
                Spro = [A'*P+P*A+opt.eps*Inx 	P*B; 
					   B'*P 						Znz];

			  	for i = 1:NMeq
			  		Spro = Spro + nu(i)*Meq(:,:,i);
			  	end
			  	for i = 1:NMi
			  		Spro = Spro + xi(i)*(alp_test^2*Mi_tilde(:,:,i) + Mi_hat(:,:,i));
			  	end
			  	Spro <= 0;
        cvx_end
        
        alp_old = alp_test;
        
        % check feasibility        
        if strcmp(cvx_status, 'Solved')
            % solve for r
            r = 1/sqrt(lam);
            if opt.verbose
                fprintf('alp_test = %2.4f: SDP feasible! r = %2.6f \n', alp_test, r);
            end
            
            % save info
            info.feasibility = true;
            info.r = r;
            info.P = P;
            info.nu = nu;
            info.xi = xi;
            info.LMI = Spro;            
            
            % alp update
            alp_min = alp_test;
            alp_test = (alp_test + alp_max)/2;
        else
            r = -Inf;
            if opt.verbose
                fprintf('alp_test = %2.4f: SDP infeasible \n', alp_test);
            end
            
            % save info
            info.feasibility = false;
            info.r = r;
            
            % alp update
            alp_max = alp_test;
            alp_test = (alp_min + alp_test)/2;
        end
        
        % update list
        list_r(cnt_iter) = r;
        list_info{cnt_iter} = info;
        
        % update iteration condition
        del_alp = alp_test - alp_old;
     end
    
    % find the largest r in the list
    [r_max, idx] = max(list_r);
    info_max = list_info{idx};    
end
            
% function info_stru = get_infoStru()
% 	info_stru = struct('feasibility',false, 'r', nan, ...
% 					'alpha', nan, 'E', nan, 'P', nan, 'nu', nan, 'xi', nan, ...
% 					'LMI', nan);
% end