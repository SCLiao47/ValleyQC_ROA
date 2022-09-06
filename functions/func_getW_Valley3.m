
% [test case1]
%   input:
%       model.nx = 2;
%       Q = diag([1,1,-1]);
%       E = eye(3);
%   output:
%       W = [2 0 0; 0 2 2; 0 2 2], [2 0 0; 0 2 -2; 0 -2 2], 
%           [2 0 2; 0 2 0; 2 0 2], [2 0 -2; 0 2 0; -2 0 2];

% [test case2]
%   input:
%       model.nx = 2;
%       Q = diag([1,1,-1]);
%       E = eye(3);
%   output:
%       W = [2 2 0; 2 2 0; 0 0 2], [2 -2 0; -2 2 0; 0 0 2], 
%           [2 0 2; 0 2 0; 2 0 2], [2 0 -2; 0 2 0; -2 0 2];


% return W of rank3 QC on nonlineaity Q
function W = func_getW_Valley3(model, Q, E)
	% system properties
	nx = model.nx;
	Znx = zeros(nx);

	% initialization
	W = repmat(Znx,[1,1,4]);

	%% Bounding procedure
	[V,D] = eig(Q);

	idx_p = find(diag(D) > 0);
	idx_n = find(diag(D) < 0);

    %% vectors nu1, nu2, nu3
    lams = zeros(3,1);
    vs = zeros(nx,3);
    if length(idx_p) == 2
        % Q has 2 positive eigenvalues and 1 negative eigenvalue
        lams(1) = D(idx_p(1), idx_p(1));
        lams(2) = D(idx_p(2), idx_p(2));
        lams(3) = D(idx_n, idx_n);
        
        vs(:,1) = V(:, idx_p(1));
        vs(:,2) = V(:, idx_p(2));
        vs(:,3) = V(:, idx_n);
    else
        % Q has 1 positive eigenvalue and 2 negative eigenvalues
        lams(1) = D(idx_n(1), idx_n(1));
        lams(2) = D(idx_n(2), idx_n(2));
        lams(3) = D(idx_p, idx_p);
        
        vs(:,1) = V(:, idx_n(1));
        vs(:,2) = V(:, idx_n(2));
        vs(:,3) = V(:, idx_p);        
        
        Q = -Q;     % invert the sign
    end
    
    nus = zeros(nx,3);
    for i = 1:3
        nus(:,i) = sqrt(abs(lams(i))) * vs(:,i);
    end
    
    %% first set of QCs, c1 = nu1 + nu3, c2 = nu1 - nu3
    c1 = nus(:,1) + nus(:,3);
    c2 = nus(:,1) - nus(:,3);
    S = nus(:,2)*nus(:,2)';
    
    W(:,:,1:2) = bounding_V3(Q, c1, c2, S, E);
    
    %% second set of QCs, c1 = nu2 + nu3, c2 = nu2 - nu3
    c1 = nus(:,2) + nus(:,3);
    c2 = nus(:,2) - nus(:,3);
    S = nus(:,1)*nus(:,1)';
    
    W(:,:,3:4) = bounding_V3(Q, c1, c2, S, E);
end

function Ws = bounding_V3(Q, c1, c2, S, E)
    invE = inv(E);
    invEsq = sqrtm(invE);
    
    beta1 = c2'*invE*c2;
    beta2 = c1'*invE*c1;
    gamma = func_lambdaMax(invEsq*(2*Q - S)*invEsq);
    
    Ws = repmat(0*E, 1, 1, 2);
    Ws(:,:,1) = beta1*(c1*c1') + gamma*S;
    Ws(:,:,2) = beta2*(c2*c2') + gamma*S;
end