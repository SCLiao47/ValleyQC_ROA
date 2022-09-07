% Scripts that implements the largest spherical ROA estimation algorithm [1] for a 
% 2-state model in [2]. This analysis includes our proposed monomial formulation and valley QCs


% [Reference]
%   [1] A. Kalur, T. Mushtaq, P. Seiler, and M. S. Hemati, “Estimating regions
% 	of attraction for transitional flows using quadratic constraints,”
% 	IEEEdfgsd
% 	Control Systems Letters, 2021.
%   [2] F. Amato, C. Cosentino, and A. Merola, “On the region of asymptotic
% 	stability of nonlinear quadratic systems,” in 2006 14th Mediterranean
% 	Conference on Control and Automation, pp. 1–5, IEEE, 2006.

clc; clear; close all;
init;

%% models
createModel = @model_2State;
createModel_mon = @model_2State_monomial;

% create model
% model = createModel();
model = createModel_mon();

%% Settings
% using spherical shape
E0 = eye(model.nx);

% [Set 1] QCs to use
CS_eQC = {};
CS_ieQC = {'CS_z'};

% [Set 2] QCs to use
CSValley_eQC = {};
CSValley_ieQC = {'CS_z';'Valley_z'};

% Optimization options
niter = 1;

options_Opt_ROA = func_getOptions_SDP_ROA(1,-2,201,-6,false);
options_Opt_ROA.verbose = true;

options_Shape.Opt_ROA = options_Opt_ROA;
options_Shape.verbose = true;


%% Run Analysis
% r_xx is the ROA estimation size

% Set 1:
disp('===[ CS ShapeIter ]===');
[r_CS, info_CS, rs_CS, infos_CS] = func_ShapeIteration(model, ...
                                                        CS_eQC, CS_ieQC, ...
                                                        E0, niter, ...
                                                        options_Shape);

% Set 2:
disp('===[ CS_valley ShapeIter ]===');
[r_CSValley, info_CSValley, rs_CSValley, infos_CSValley] = func_ShapeIteration(model, ...
                                                        CSValley_eQC, CSValley_ieQC, ...
                                                        E0, niter, ...
                                                        options_Shape);                                             

%% save data
save(['Data\ROA_2State\', datestr(now, 'yyyymmdd_HHMM')]);

%% plotting
plotting_2State;

%% Result summary
fprintf('[Analysis Results] \n');
fprintf('Set 1: r = %.4f \n', r_CS);
fprintf('Set 2: r = %.4f \n', r_CSValley);