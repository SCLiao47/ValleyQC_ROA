


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

%% Analysis cases
% Quadratic constraints
numCases = 5;
ieQCs = cell(numCases,1);

ieQCs{1} = {'CS_z'; 'CS_phi'};
ieQCs{2} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'};
ieQCs{3} = {'CS_z'; 'CS_phi'; 'Valley3_phi'};
ieQCs{4} = {'CS_z'; 'CS_phi'; 'CrossP_z'; 'CrossP_neg'};
ieQCs{5} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'; 'Valley3_phi'; 'CrossP_z'; 'CrossP_neg'};

QCcases.ieQC = ieQCs;

%%
load('Random3States\Exps2_20220320_0025_EECS.mat');

model = ExpList(2).model;
model.Mi.CrossP_neg = @(E) func_getQCs(model,'Mi_CrossProduct_z',E);

exp = run_analysis(model, QCcases, opt);

%%
metrics = get_metrics(exp.r)