% try to find [numExp] of model that all other sets of QC perform better
% than the baseline QC (CSQC)

clc; clear; close all;
init;

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

opt.E0 = E0;
opt.niter = 1;
opt.options_Shape = options_Shape;

%% simulation cases
% Quadratic constraints
numCases = 5;
ieQCs = cell(numCases,1);

ieQCs{1} = {'CS_z'; 'CS_phi'};
ieQCs{2} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'};
ieQCs{3} = {'CS_z'; 'CS_phi'; 'Valley3_phi'};
ieQCs{4} = {'CS_z'; 'CS_phi'; 'CrossP_z'};
ieQCs{5} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'; 'Valley3_phi'; 'CrossP_z'};

QCcases.ieQC = ieQCs;

%% Recording setupt
numExp = 5;
countExp = 0;

ExpList = [];

%% start simulation
modelCount = 0;

while countExp < numExp
    modelCount = modelCount + 1;
    disp(['--- Model #', num2str(modelCount)]);
    
    model = createModel();
    
    % run analysis
    exp = run_analysis(model, QCcases, opt);

    % get metrics and store analysis
    metrics = get_metrics(exp.r);
    exp.metrics = metrics;
    
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

if ~exist('Data\Random3States')
    mkdir('Data\Random3States')
end
save(['Data\Random3States\Exps',num2str(numExp), '_', datestr(now, 'yyyymmdd_HHMM')]);
    