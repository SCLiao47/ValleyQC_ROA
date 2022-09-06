function Mi = func_getCrossQC_neg(model, E)
    % system properties
    nx = model.nx;
    nz = model.nz;
    
    nxz = nx + nz;
    
    Qis = model.Qis;
    
    %% check which pairs to generate QC
    w_pairs = [];
    
    for p = 1:nz
        for q = p+1:nz
            % check power of x of Qi
            cp = get_XpowerInQ(Qis(:,:,p));
            
            % check power of x of Qj
            cq = get_XpowerInQ(Qis(:,:,q));
            
            %% [check cases]
            QCflag = false;
            
            idx_xExponent = find(cp+cq);
            
            switch length(idx_xExponent)
                case 1      % case xi^4
                    warning(['Something went wrong in Qs! Check w: ',...
                                    num2str(p),', ', num2str(q)]);
                    
                case 2      % case xi^2*xj^2 or xi^3*xj
                    if length(find(cp+cq == 2)) == 2    % case xi^2*xj^2
                        QCflag = false;
                    else                                % case xi^3*xj
                        new_pair = struct('p',p, 'q',q);
                        
                        idx_power3 = find(cp+cq == 3);
                        idx_power1 = find(cp+cq == 1);
                        
                        new_pair.i = idx_power3;
                        new_pair.j = idx_power3;
                        new_pair.k = idx_power1;
                        
                        QCflag = true;
                    end
                    
                case 3      % case xi^2*xj*xk
                    new_pair = struct('p',p, 'q',q);
                    
                    new_pair.i = find(cp+cq == 2);
                    
                    idx_power1 = find(cp+cq == 1);
                    new_pair.j = idx_power1(1);
                    new_pair.k = idx_power1(2);
                    
                    QCflag = true;
                
                case 4      % case xi*xj*xk*xl: skip
                    QCflag = false;
                    
                otherwise
                    warning(['Something went wrong in Qs! Check w: ',...
                                    num2str(p),', ', num2str(q)]);
            end           
            
            % add into list
            if QCflag
                w_pairs = cat(1, w_pairs, new_pair);
            end
        end
    end
    
    %% generate QCs for each pair in w_paris
    numQC = length(w_pairs);
    
    % initialize QCs
    Mi.tilde = zeros(nxz, nxz, 2*numQC);
    Mi.hat = zeros(nxz, nxz, 2*numQC);
    
    for idx = 1:numQC
        pair = w_pairs(idx);
        
        M = func_getCrossProductQC_negative(model,E, pair.p, pair.q, ...
                                    pair.i, pair.j, pair.k);
                                
        idx_M = 2*(idx-1) + [1,2];
        Mi.tilde(:,:,idx_M) = M.tilde;
        Mi.hat(:,:,idx_M) = M.hat;
    end
end


function c = get_XpowerInQ(Q)
    c = sum(Q,1)' + sum(Q,2); 
end