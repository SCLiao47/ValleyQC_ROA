% This function solves the SDP (8) in the paper for ROA estimate given set of QCs within local region x'Ex. 
% Note that the local region is predefined outside of this function.



function [r_max, info_max, grid_r, grid_info] = SDPOpt_ROAEstimation(model, E, Meq, Mi, options)
	% check input argument
	if nargin < 5
        options = func_getOptions_SDP_ROA();
	end

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

    %% solve SDP for over a grid of alpha
    % initialize grid of alpha 
%     grid_alp = logspace(log10(options.alp_min), log10(options.alp_max), options.Nalp);
    grid_alp = [logspace(log10(options.alp_min), log10(options.alp_max), options.Nalp-1) 1];

    % initialize grid of r and info w.r.t. alpha values
    grid_r = -Inf(options.Nalp, 1);

    grid_info = cell(options.Nalp, 1);
    [grid_info{:}] = deal(get_infoStru());

    % initialize info variable for largest r*
    r_max = -Inf;
    info_max = get_infoStru();

    for idx_alp = 1:options.Nalp
    	alp = grid_alp(idx_alp);
        
        % initialize info
        info = get_infoStru();
        info.alpha = alp;
        info.E = E;

    	% [convex SDP problem]
    	% lam = 1/r^2;
		cvx_begin sdp quiet
		cvx_solver mosek
			variable P(nx,nx)  		semidefinite;
			variable nu(NMeq,1);
			variable xi(NMi,1)      nonnegative;
			variable lam 	 		nonnegative;

			minimize(lam);

			subject to
				% set containment
				P <= lam*Inx;
				E/alp^2 <= P;

				% Lyapunov stability condition
				Spro = [A'*P+P*A+options.eps*Inx 	P*B; 
					   B'*P 						Znz];

			  	for i = 1:NMeq
			  		Spro = Spro + nu(i)*Meq(:,:,i);
			  	end
			  	for i = 1:NMi
			  		Spro = Spro + xi(i)*(alp^2*Mi_tilde(:,:,i) + Mi_hat(:,:,i));
			  	end
			  	Spro <= 0;
	  	cvx_end

	  	% check if feasible
	  	if strcmp(cvx_status, 'Solved')
	  		isFesible = true;

	  		r = 1/sqrt(lam);

	  		if options.verbose
	  			fprintf('alpha = %2.4f: SDP feasible! r = %2.6f \n', alp, r);
	  		end

	  		% setting info

            
	  		info.feasibility = true;
	  		info.r = r;
	  		info.P = P;
	  		info.nu = nu;
	  		info.xi = xi;
            info.LMI = Spro;
	  	else
	  		r = -Inf;

	  		if options.verbose
                fprintf('alpha = %2.4f: SDP infeasible \n', alp);
		  	end
	  	end

	  	% storing info grid
	  	grid_r(idx_alp) = r;
	  	grid_info{idx_alp} = info; 

	  	% compare to r_max
	  	if r > r_max
	  		r_max = r;
	  		info_max = info;
	  	end
  	end
end

function Spro = get_SProcedure(model, E, alp, Meq, Mi, NMeq, NMi, P, nu, xi)
    % system dynamics
	A = model.A;
	B = model.B;
	nx = model.nx;
	nz = model.nz;
    
    % S-procedure
    Spro = [A'*P+P*A+options.eps*Inx 	P*B; 
					   B'*P             Znz];

    for i = 1:NMeq
        Spro = Spro + nu(i)*Meq(:,:,i);
    end
    for i = 1:NMi
        Spro = Spro + xi(i)*(alp^2*Mi.tilde(:,:,i) + Mi.hat(:,:,i));
    end 
end


function info_stru = get_infoStru()
	info_stru = struct('feasibility',false, 'r', nan, ...
					'alpha', nan, 'E', nan, 'P', nan, 'nu', nan, 'xi', nan, ...
					'LMI', nan);
end