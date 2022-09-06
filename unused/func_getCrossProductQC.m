

function M = func_getCrossProductQC(model, E, p, q, i, j, k)
    % system dimension
    nx = model.nx;
    Inx = eye(nx);
    Znx = zeros(nx);
    
    nz = model.nz;
    Inz = eye(nz);
    Znz = zeros(nz);
    
    nxz = nx+nz;
    
    % forming QC sturctures
    M.tilde = zeros(nxz,nxz,4);
    M.hat = zeros(nxz,nxz,4);
    
    % auxillary vector
    ei = get_StandardBasis(nx, i);
    ej = get_StandardBasis(nx, j);
    ek = get_StandardBasis(nx, k);
    
    cp = ej + ek;
    cn = ej - ek;
    
    ebarp = get_StandardBasis(nz, p);
    ebarq = get_StandardBasis(nz, q);
    
    % forming matrices
    invE = inv(E);
    
    Spq = ebarp*ebarq' + ebarq*ebarp';
    Wp1 = (ei'/E*ei)*(cp*cp');
    Wp2 = (cp'/E*cp)*(ei*ei');
    Wn1 = (ei'/E*ei)*(cn*cn');
    Wn2 = (cn'/E*cn)*(ei*ei');
    
    % setting output
    M.tilde(:,:,1) = blkdiag(Wp1, Znz);
    M.tilde(:,:,2) = blkdiag(Wp2, Znz);
    M.tilde(:,:,3) = blkdiag(Wn1, Znz);
    M.tilde(:,:,4) = blkdiag(Wn2, Znz);
    
    M.hat(:,:,1) = blkdiag(Znx, -Spq);
    M.hat(:,:,2) = blkdiag(Znx, -Spq);
    M.hat(:,:,3) = blkdiag(Znx, Spq);
    M.hat(:,:,4) = blkdiag(Znx, Spq);
end
    
function ei = get_StandardBasis(n, i)
    ei = zeros(n,1);
    ei(i) = 1;
end