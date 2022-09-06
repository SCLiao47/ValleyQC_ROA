% Function return the option argument for SDPOpt_ROAEstimation

function options = func_getOptions_SDP_ROA(alp_max, alp_min, Nalp, eps, verbose)
    if nargin == 0
        alp_max = 0;
        alp_min = -2;
        Nalp = 20;
        eps = -6;
        verbose = false;
    end
    
    options.alp_max = power(10,alp_max);
    options.alp_min = power(10,alp_min);
    options.Nalp = Nalp;

    options.eps = power(10,eps);

    options.verbose = verbose;
end