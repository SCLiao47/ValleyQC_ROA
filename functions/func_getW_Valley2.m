

% return W of Valley QC on nonlineaity Q 
function [W, L] = func_getW_Valley2(model, Q, E)
	% ### TODOs
	%  	- check sign-indefinite
    eig_eps = 1e-6;

	% system properties
	nx = model.nx;
	Znx = zeros(nx);

	% initialization
	W = repmat(Znx,[1,1,2]);

	%% Bounding procedure
	[V,D] = eig(Q);

	idx_p = find(diag(D) > eig_eps); 
	idx_n = find(diag(D) < -eig_eps);

    if or(length(idx_p) > 1, length(idx_n) > 1)
        warning('Q is not rank 2 sign-indefinite. Take the last idx_p and first idx_n.');
        
        idx_p = idx_p(end);
        idx_n = idx_n(1);
    end
    
	lamb_p = D(idx_p,idx_p);
	lamb_n = D(idx_n,idx_n);
	v_p = V(:,idx_p);
	v_n = V(:,idx_n);

	nu_p = sqrt(lamb_p)*v_p;
	nu_n = sqrt(-lamb_n)*v_n;

    c1 = nu_p + nu_n;
    c2 = nu_p - nu_n;

    W1 = (c2'*inv(E)*c2) * c1*c1';
    W2 = (c1'*inv(E)*c1) * c2*c2';
    W = cat(3,W1,W2);
    
    % for new xi-E algorithm
    L = cat(3, c1*c2', c2*c1');
end