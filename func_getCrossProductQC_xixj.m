

function M = func_getCrossProductQC_xixj(model, E, p, q, i, j)
    % system dimension
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % forming QC sturctures
    M.tilde = zeros(nxz,nxz,2);
    M.hat = zeros(nxz,nxz,2);
    
    % auxillary vector
    ei = get_StandardBasis(nx, i);
    ej = get_StandardBasis(nx, j);
    
    ebarp = get_StandardBasis(nz, p);
    ebarq = get_StandardBasis(nz, q);
    
    % forming matrices
    invE = inv(E);
    
    Spq = ebarp*ebarq' + ebarq*ebarp';
    W1 = (ei'/E*ei)*(ej*ej');
    W2 = (ej'/E*ej)*(ei*ei');
    
    % setting output
    M.tilde(:,:,1) = blkdiag(W1, Znz);
    M.tilde(:,:,2) = blkdiag(W2, Znz);
    
    M.hat(:,:,1) = blkdiag(Znx, -Spq);
    M.hat(:,:,2) = blkdiag(Znx, -Spq);
end
    
function ei = get_StandardBasis(n, i)
    ei = zeros(n,1);
    ei(i) = 1;
end