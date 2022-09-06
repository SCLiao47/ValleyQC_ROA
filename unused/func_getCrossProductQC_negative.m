

function M = func_getCrossProductQC_negative(model, E, p, q, i, j, k)
    % system dimension
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % auxillary vector
    ei = get_StandardBasis(nx, i);
    ej = get_StandardBasis(nx, j);
    ek = get_StandardBasis(nx, k);
    
    c = ej - ek;
    
    ebarp = get_StandardBasis(nz, p);
    ebarq = get_StandardBasis(nz, q);
    
    % forming matrices
    invE = inv(E);
    
    Spq = ebarp*ebarq' + ebarq*ebarp';
    W1 = (ei'/E*ei)*(c*c');
    W2 = (c'/E*c)*(ei*ei');
    
    % forming QC sturctures
    M.tilde = zeros(nxz,nxz,2);
    M.hat = zeros(nxz,nxz,2);
    
    % setting output
    M.tilde(:,:,1) = blkdiag(W1, Znz);
    M.tilde(:,:,2) = blkdiag(W2, Znz);
    
    M.hat = repmat(blkdiag(Znx, Spq), 1,1,2);
end
    
function ei = get_StandardBasis(n, i)
    ei = zeros(n,1);
    ei(i) = 1;
end