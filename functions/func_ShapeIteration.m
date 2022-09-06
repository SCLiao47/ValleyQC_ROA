% Function to run shaper iteration algorithm (Algorithm A) in [~]



function [r_final, info_final, r_list, info_list] = func_ShapeIteration(model, eQC_names, ieQC_names, E0, niter, options)
    % check options
    if nargin < 6
        options.Opt_ROA = func_getOptions_SDP_ROA();
        options.verbose = false;
    end

    % initialize output
    r_list = zeros(niter,1);
    info_list = cell(niter,1);

    % initialize shape
    E = E0;

    for iter = 1:niter
        if options.verbose
%             disp('=======================');
            fprintf('--- Shaper iteration #%2i \n',iter);
        end
        
%         if iter > 1
%             options.Opt_ROA.alp_min = 1;
%             
%             options.Opt_ROA.alp_max = 1 + power(2, -iter);
%         end
        
        % get M0 and Mi
        Meq = func_stackingQCs(model.Meq, eQC_names, 'eQC',E);
        Mi = func_stackingQCs(model.Mi, ieQC_names, 'ieQC',E);        
        
        % run SDP
        [ri, info_i] = SDPOpt_ROAEstimation(model, E, ...
                                            Meq, Mi, ...
                                            options.Opt_ROA);
                                        
        if options.verbose
            fprintf('    Niter = %2i: ri = %2.4f \n', iter, ri);
        end
                                        
        % store output
        r_list(iter) = ri;
        info_list{iter} = info_i;
                                        
        % update E
        E = info_i.P;
    end
    
    % set outputs
    r_final = r_list(niter);
    info_final = info_list{niter};
end