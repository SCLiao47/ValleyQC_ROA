%% axuliary function
function metrics = get_metrics(rCases)
    % options
%     impRate = 1e-3;
    impRate = 0.05;

    % compute metrics
    numCases = length(rCases);

    base = rCases(1);
    rate = (rCases - base) / base;
        
    % maximum improvement
    metrics.maxImp = max(rate);
    
    % average improvement
    metrics.meanImp = sum(rate)/(numCases-1);
    
    % number of improved cases
    isImp = rate >= impRate;
    metrics.numImp = length(find(isImp));
    
    % if all QC makes improve?
    isAllImp = all(rCases(end) > rCases(1:end-1));
    metrics.isAllImp = isAllImp;
end