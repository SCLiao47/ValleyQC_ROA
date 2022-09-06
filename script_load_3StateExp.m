
% loading model
load('Data\Random3States\Exps8_20220321_2235_EECS.mat')
exp = ExpList(4);
model = exp.model;

% setting QCs to analyze
ieQCs{1} = {'CS_z'; 'CS_phi'};
ieQCs{2} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'};
ieQCs{3} = {'CS_z'; 'CS_phi'; 'Valley3_phi'};
ieQCs{4} = {'CS_z'; 'CS_phi'; 'CrossP_z'};
ieQCs{5} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'; 'Valley3_phi'};
ieQCs{6} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'; 'CrossP_z'};
ieQCs{7} = {'CS_z'; 'CS_phi'; 'Valley3_phi'; 'CrossP_z'};
ieQCs{8} = {'CS_z'; 'CS_phi'; 'Valley2_z'; 'Valley2_phi'; 'Valley3_phi'; 'CrossP_z'};

QCcases.ieQC = ieQCs;

% setting options
opt.options_Shape.Opt_ROA.Nalp = 200;

% run analysis
new_exp = run_analysis(model, QCcases, opt);

% display result
new_exp.r
