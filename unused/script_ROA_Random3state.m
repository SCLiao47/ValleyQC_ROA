


% clc; clear; close all;

%% model 
createModel = @model_Random3States;

%% options

% ROA estimation
options_OptROA = func_getOptions_SDP_ROA(1,-3,50,-6,false);

% shape iteration
E0 = eye(3);
niter = 1;
options_Shape.Opt_ROA = options_OptROA;
options_Shape.verbose = false;


%% simulation cases
% Quadratic constraints
numCases = 4;
ieQCs = cell(numCases,1);

ieQCs{1} = {'CS_z'; 'CS_phi'};
ieQCs{2} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'};
ieQCs{3} = {'CS_z'; 'CS_phi'; 'Valley3_phi'};
ieQCs{4} = {'CS_z'; 'CS_phi'; 'CrossP_z'};
% ieQCs{5} = {'CS_z'; 'CS_phi'; 'Valley3_phi'; 'CrossP_z'};

%% Recording setupt

numExp = 1;
countExp = 0;

ExpList = [];

%% start simulation

modelCount = 0;

while countExp < numExp
    modelCount = modelCount + 1;
    disp(['--- Model #', num2str(modelCount)]);
    
    model = createModel();

    rCases = zeros(numCases, 1);
    infoCases = cell(numCases, 1);

    for idx_case = 1:numCases
        disp(['  - QC case #', num2str(idx_case)]);

        [rCases(idx_case), infoCases{idx_case}, ~, ~] = func_ShapeIteration(model, ...
                                                    [], ieQCs{idx_case}, ...
                                                    E0, niter, ...
                                                    options_Shape);
    end     
    
    % store exp
    exp.model = model;
    exp.r = rCases;
    exp.info = infoCases;

    % get metrics and store analysis
    metrics = get_metrics(rCases);
    
    disp([' - numImp: ', num2str(metrics.numImp)]);
    if metrics.numImp == numCases-1
        ExpList = [ExpList; exp];
        
        disp(' - add to the list!');
        disp([' - maxImp: ', num2str(metrics.maxImp)]);
        
        countExp = countExp + 1;
    end
    
    fprintf('\n');
end

%%
disp('Model that all QC cases perform better than baseline(CSQC) found!');
disp('Saving the models to [~/Random3States/]');

if ~exist('Random3States')
    mkdir('Random3States')
end
save(['Random3States\Exps',num2str(numExp), '_', datestr(now, 'yyyymmdd_HHMM')]);

%% axuliary function
function metrics = get_metrics(rCases)
    % options
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
end
    