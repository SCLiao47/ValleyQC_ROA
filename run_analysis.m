


function exp = run_analysis(model, QCcases, opt)
    ieQCs = QCcases.ieQC;

    numCases = size(ieQCs,1);
    
    rCases = zeros(numCases, 1);
    infoCases = cell(numCases, 1);

    for idx_case = 1:numCases
        disp(['  - QC case #', num2str(idx_case)]);

        [rCases(idx_case), infoCases{idx_case}, ~, ~] = func_ShapeIteration(model, ...
                                                    [], ieQCs{idx_case}, ...
                                                    opt.E0, opt.niter, ...
                                                    opt.options_Shape);
    end     
    
    % store exp
    exp = [];
    exp.model = model;
    exp.r = rCases;
    exp.info = infoCases;
end