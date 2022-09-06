

function [r, info] = func_ROAnewAlgorithm(model, eQC_names, ieQC_names, E0, opt)
    % check potions
    if ~isfield(opt, 'r_tol')
        opt.r_tol = 0.01;
    end
    if ~isfield(opt, 'Opt_ROA')
        opt.Opt_ROA = func_getOptions_SDP_ROA(1,-2,201,-6,false);
    end
    if ~isfield(opt, 'alp_bisection')
        opt.alp_bisection = true;
    end
    
%     opt.opt_ROA.verbose = true;
    
    %% [initial feasible solution] E = E0
    % get Meq and Mi
    Meq = func_stackingQCs(model.Meq, eQC_names, 'eQC',E0);
    Mi = func_stackingQCs(model.Mi, ieQC_names, 'ieQC',E0); 
    
    % run SDP
    if opt.alp_bisection
        [r_1, info_1, num_iter_alpha] = SDPOpt_ROA_ALP_bisection(model, E0, Meq, Mi, opt.Opt_ROA);
    else
        [r_1, info_1] = SDPOpt_ROAEstimation(model, E0, Meq, Mi, opt.Opt_ROA);
        num_iter_alpha = 200;
    end
        
    % terminate function if infeasible
    if opt.verbose
        if ~info_1.feasibility
            fprintf('Initialization with E=I failed. \n');
            return;
        else
            fprintf('[Initialization with E=I] r = %2.4f with alpha = %2.4f \n', r_1, info_1.alpha);
        end
        
        fprintf('    # iteration: %i \n', num_iter_alpha);
    end
    
    %% Alternative solving different variables
    fprintf('------[ E-xi iteration ]--- \n');
    
    % initialization
    info = info_1;    
    r = r_1;
    del_r = inf;
    num_iter_Exi = 0;
    
    while abs(del_r) > opt.r_tol
        num_iter_Exi = num_iter_Exi + 1;
        
        % solving for {P,r,Etild}
        [r_E, info_E] = SDPOpt_ROA_Etile(model, info.xi, Meq, Mi, opt.Opt_ROA);
 
        
        % solving for {P,r,xi}
        Meq = func_stackingQCs(model.Meq, eQC_names, 'eQC',info_E.Etilde);
        Mi = func_stackingQCs(model.Mi, ieQC_names, 'ieQC',info_E.Etilde); 
        
        [r_xi, info_xi] = SDPOpt_ROA_xi(model, info_E.Etilde, Meq, Mi, opt.Opt_ROA);        
        
        % checking resutls
        if r > r_E + 1e-8
            warning('E-step is not accurate!');
        end
        if r_E > r_xi + 1e-8
            warning('Xi-step is not accurate!');
        end
                
        % update loop condition
        info = info_xi;
        del_r = (r_xi - r)/r;
        
        r = r_xi;
    end 
    if opt.verbose
        fprintf('[E-xi iteration] r = %2.6f \n', r_xi);
        fprintf('    # iteration: %i \n', num_iter_Exi);
    end
end