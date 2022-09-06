% The function returns the specified type of QC of the given model in local region x'Ex 


% [QC_type]


% [Note]
% 	Mi(alpha) = M_hat + alpha^2*M_tilde

function M = func_getQCs(model, QC_type, E)
	% ### TODOs
	% 	- check E is positive definite

	% dimension of the model
	nx = model.nx;
	nz = model.nz;

	nxz = nx+nz;

    % utility matrics
%     Inx = eye(nx);
    Znx = zeros(nx);
%     Znx1 = zeros(nx,1);
    
    Inz = eye(nz);
    Znz = zeros(nz);
%     Znz1 = zeros(nz,1);
    
%     Znxnz = zeros(nx,nz);

    % dynamics of the model
    Qis = model.Qis;
    evalues = zeros(nx,nz);	 % i-th row is the eigenvalues of i-th Qi	
	eig_max = zeros(nz,1);
	eig_min = zeros(nz,1);

	for i = 1:nz
		Qi = Qis(:,:,i);

		evalues(:,i) = eig(Qi);
		eig_max(i) = max(evalues(:,i));
		eig_min(i) = min(evalues(:,i));
	end

	%% return M of specified type of QC
	switch QC_type
		%% [Meq]
		case 'Meq_lossless'
			% x*B*z = 0
			M = [Znx model.B/2; model.B'/2 Znz];

		%% [Mi]
		case 'Mi_CS_z' 		% Cauchy-Schwarz for each nonlinearity zi
            % zi^2 = (x'*Qi*x)^2
                 % = (w'*Qi_hat*w)^2,				(w=E^(0.5)x)
            %      <= ||w||^2 * ||Qi_hat*w||^2
            %      <= alpha^2 * x'*(Qi*E^-1*Qi)*x
            %      = x' * (alpha^2*Qi*E^-1*Qi) * x

            % initialization
            M.tilde = blkdiag(Znx,Znz);
            M.hat = blkdiag(Znx,Znz);
            M.L = Znx;

            % for each zi = x'*Qis(:,:,i)*x
            for i = 1:nz
                Qi = Qis(:,:,i);
            	ei = Inz(:,i);

            	M.tilde(:,:,i) = blkdiag(Qi/E*Qi, Znz);
            	M.hat(:,:,i) = blkdiag(Znx, -ei*ei');
                
                M.L(:,:,i) = Qi;
            end

        case 'Mi_CS_phi'	% Cauchy-Schwarz for each nonlinearity in phi = B*z;
            if ~isfield(model, 'Qhats')
            	error('QC [Mi_CS_phi] is not implemented for this model yet!');
            end
            
            % initialization
            M.tilde = blkdiag(Znx,Znz);
            M.hat = blkdiag(Znx,Znz);
            M.L = Znx;
            
            % for each B(i,:)*w
            for i = 1:nx
                Qhat = model.Qhats(:,:,i);
                b = model.B(i,:)';
                
                M.tilde(:,:,i) = blkdiag(Qhat/E*Qhat, Znz);
                M.hat(:,:,i) = blkdiag(Znx, -b*b');
                
                M.L(:,:,i) = Qhat;
            end
                

        case 'Mi_Valley_z' 	% Valley QCs for each nonlinearity zi
        	% ### TODOs
        	% 	- check only 2 non-zero eigenvalues

        	% check which zi to bound (sign-indefinite Qi)
        	ids_indefinite = find(eig_max.*eig_min < 0);
        	n_indefinite = size(ids_indefinite,1);

        	% initialize  	(each sign-indefinite zi has 2 valley QC)
        	M.tilde = zeros(nxz,nxz, n_indefinite*2);
        	M.hat = zeros(nxz,nxz, n_indefinite*2);
            M.L = zeros(nx,nx, n_indefinite*2);

        	% for each sign-indefinite zi = x'*Qis(:,:,i)*x
            for i = 1:n_indefinite
                k = ids_indefinite(i); 	 	% k-th element of z
                ek = Inz(:,k);

                % call bounding function to get W
                [W, L] = func_getW_Valley2(model, Qis(:,:,k), E);

                % assign outputs
                for j = 1:2
                    M.tilde(:,:,(i-1)*2+j) = blkdiag(W(:,:,j), Znz);
                    M.hat(:,:,(i-1)*2+j) = blkdiag(Znx, -ek*ek');
                    
                    M.L(:,:,(i-1)*2+j) = L(:,:,j);
                end
            end
            
        case 'Mi_Valley_phi'    % Valley QC for each nonlinearity bi'*w
            ids_Valley2_phi = [];
            
            % check how many rank-2 sign-indefinite Qhat
            for i = 1:nx
                Qhat = model.Qhats(:,:,i);
                
                % check rank 2 to initiate the process
                info_Qhat = get_matrixInfo(Qhat);
                if strcmp(info_Qhat.class, 'rank2_indefinite')
                    ids_Valley2_phi = [ids_Valley2_phi; i];
                end
            end
               
            numQC = length(ids_Valley2_phi);
            
            % initialize
            M.tilde = zeros(nxz,nxz, numQC*2);
        	M.hat = zeros(nxz,nxz, numQC*2);
            M.L = zeros(nx, nx, numQC*2);

            % assign Valley QC for each rank-2, sign-indefinte Qhat
            for i = 1:numQC
                idx = ids_Valley2_phi(i);
                b = model.B(idx,:)';
                
                % call bounding function to get W
                [W, L] = func_getW_Valley2(model, model.Qhats(:,:,idx), E);
                
                % assign outputs
                for j = 1:2
                    idx_M = (i-1)*2+j;
                    M.tilde(:,:,idx_M) = blkdiag(W(:,:,j), Znz);
                    M.hat(:,:,idx_M) = blkdiag(Znx, -b*b');
                    
                    M.L(:,:,idx_M) = L(:,:,j);
                end
            end
            
        case 'Mi_Valley3_phi'
            ids_Valley3_phi = [];
            
            % check how many rank-3, sign-indefinite Qhat
            for i = 1:nx
                Qhat = model.Qhats(:,:,i);
                
                info_Qhat = get_matrixInfo(Qhat);
                if strcmp(info_Qhat.class, 'rank3_indefinite')
                    ids_Valley3_phi = [ids_Valley3_phi; i];
                end
            end
            
            numQC = length(ids_Valley3_phi);
            
            % initialize    (each Qhat of the class has 4 rank-3 QC)
            M.tilde = zeros(nxz,nxz, numQC*4);
            M.hat = zeros(nxz,nxz, numQC*4);
 
            for i = 1:numQC
                idx = ids_Valley3_phi(i);
                b = model.B(idx,:)';
                
                % call bounding function to get W
                W = func_getW_Valley3(model, model.Qhats(:,:,idx), E);
                
                % assign outputs
                for j = 1:4
                    idx_M = (i-1)*4+j;
                    M.tilde(:,:,idx_M) = blkdiag(W(:,:,j), Znz);
                    M.hat(:,:,idx_M) = blkdiag(Znx, -b*b');
                end
            end
            
        case 'Mi_CrossProduct_z'
            M = func_getCrossQC_z(model,E);
            
        case 'Mi_CrossProduct_neg'
            error('Mi_CrossProduct_neg is now included intp Mi_CrossProduct_z. Please use Mi_CrossProduct_z only');
%             M = func_getCrossQC_neg(model,E);
            
		otherwise
            error(['QC type [',QC_type, '] not defined!']);
	end
end

function info = get_matrixInfo(Q)
    % options
    eig_eps = 1e-6;

    % check rank
    info.rank = rank(Q);

    % check sign-definiteness
    eValues = eig(Q);
    
    idx_p = find(eValues > eig_eps);
    idx_n = find(eValues < -eig_eps);
    
    info.eig_np = length(idx_p);
    info.eig_nn = length(idx_n);
    
    if and(info.eig_np >= 1, info.eig_nn >= 1)
        info.sign = false;
    else
        info.sign = true;
    end
    
    % classify the matrix
    info.class = 'none';
    if and(info.rank == 2, ~info.sign)
        info.class = 'rank2_indefinite';
    end
    
    if and(info.rank == 3, ~info.sign)
        info.class = 'rank3_indefinite';
    end
end    
    