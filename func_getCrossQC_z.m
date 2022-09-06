function Mi = func_getCrossQC_z(model, E)
    % system properties
    nx = model.nx;
    nz = model.nz;
    
    nxz = nx + nz;
    
    Qis = model.Qis;
    
    %% check which pairs to generate QC
    w_pairs_ijk = [];       % xi^2*xj*xk
    w_pairs_ijsq = [];      % xi^2*xj^2
    
    for p = 1:nz
        for q = p+1:nz
            % check power of x of Qi
            cp = get_XpowerInQ(Qis(:,:,p));
            
            % check power of x of Qj
            cq = get_XpowerInQ(Qis(:,:,q));
            
            %% [check cases]
            QCflag_ijk = false;         % xi^2*xj*xk
            QCflag_ijsq = false;        % xi^2*xj^2
            
            idx_xExponent = find(cp+cq);
            
            switch length(idx_xExponent)
                case 1      % case xi^4
                    warning(['Something went wrong in Qs! Check w: ',...
                                    num2str(p),', ', num2str(q)]);
                    
                case 2      % case xi^2*xj^2 or xi^3*xj
                    if length(find(cp+cq == 2)) == 2    % case xi^2*xj^2
                        QCflag_ijsq = true;
                    else                                % case xi^3*xj
                        new_pair = struct('p',p, 'q',q);
                        
                        idx_power3 = find(cp+cq == 3);
                        idx_power1 = find(cp+cq == 1);
                        
                        new_pair.i = idx_power3;
                        new_pair.j = idx_power3;
                        new_pair.k = idx_power1;
                        
                        QCflag_ijk = true;
                    end
                    
                case 3      % case xi^2*xj*xk
                    new_pair = struct('p',p, 'q',q);
                    
                    new_pair.i = find(cp+cq == 2);
                    
                    idx_power1 = find(cp+cq == 1);
                    new_pair.j = idx_power1(1);
                    new_pair.k = idx_power1(2);
                    
                    QCflag_ijk = true;
                
                case 4      % case xi*xj*xk*xl: skip
                    QCflag_ijk = false;
                    QCflag_ijsq = false;
                    
                otherwise
                    warning(['Something went wrong in Qs! Check w: ',...
                                    num2str(p),', ', num2str(q)]);
            end           
            
            % add into list
            if QCflag_ijk
                w_pairs_ijk = cat(1, w_pairs_ijk, new_pair);
            end
            
            if QCflag_ijsq
                w_pairs_ijsq = cat(1, w_pairs_ijsq, new_pair);
            end
        end
    end
    
    %% generate QCs for each pair in w_paris
    numQC_ijk = length(w_pairs_ijk);
    numQC_ijsq = length(w_pairs_ijsq);
    
    % initialize QCs
    Mi.tilde = zeros(nxz, nxz, 4*numQC_ijk + 2*numQC_ijsq);
    Mi.hat = zeros(nxz, nxz, 4*numQC_ijk + 2*numQC_ijsq);
    
    for idx = 1:numQC_ijk
        pair = w_pairs_ijk(idx);
        
        M = func_getCrossProductQC(model,E, pair.p, pair.q, ...
                                    pair.i, pair.j, pair.k);
                                
        idx_M = 4*(idx-1) + [1,2,3,4];
        Mi.tilde(:,:,idx_M) = M.tilde;
        Mi.hat(:,:,idx_M) = M.hat;
    end
    
    for idx = 1:numQC_ijsq
        pair = w_pairs_ijsq(idx);
        
        M = func_getCrossProductQC_xixj(model, E, pair.p, pair.q, ...
                                        pair.i, pair.j);
                                    
        idx_M = 2*(idx-1) + [1,2];
        Mi.tilde(:,:,idx_M) = M.tilde;
        Mi.hat(:,:,idx_M) = M.hat;      
    end
end


function c = get_XpowerInQ(Q)
    c = sum(Q,1)' + sum(Q,2); 
end