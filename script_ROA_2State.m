% Scripts that implements the largest spherical ROA estimation algorithm [1] for a 
% 2-state model in [2]. This analysis includes our proposed monomial formulation and valley QCs


% [Reference]
%   [1] A. Kalur, T. Mushtaq, P. Seiler, and M. S. Hemati, “Estimating regions
% 	of attraction for transitional flows using quadratic constraints,” IEEE
% 	Control Systems Letters, 2021.
%   [2] F. Amato, C. Cosentino, and A. Merola, “On the region of asymptotic
% 	stability of nonlinear quadratic systems,” in 2006 14th Mediterranean
% 	Conference on Control and Automation, pp. 1–5, IEEE, 2006.

% clc; clear; close all;

%% models and settins
createModel = @model_2State;
createModel_mon = @model_2State_monomial;

%% Largest spherical ROA estimation
options_Opt_ROA = func_getOptions_SDP_ROA(1,-2,201,-6,false);
% options_Opt_ROA = func_getOptions_SDP_ROA(1,-2,20,-6,false);
% options_Opt_ROA = func_getOptions_SDP_ROA(log10(2.9), log10(2.74), 200, -6, false);

% options_Opt_ROA.Nalp = 1;
options_Opt_ROA.verbose = true;

% create model
% model = createModel();
model = createModel_mon();

% using spherical shape
E0 = eye(model.nx);

%% Checking implementation to match Liu & Gayme

% E = E0;
% 
% % [Original: CSQC]
% disp('===[ CS ]===');
% [R0max_CS, info_CS, R0s_CS] = SDPOpt_ROAEstimation(model, E,...
%                                         [], model.Mi.CS_z(E), ...
%                                         options_Opt_ROA);
% 
% disp('===[ CS Bisection ]===');
% [R0max_CS_bi, info_CS_bi, R0s_CS_bi] = SDPOpt_ROAEst_bisection(model, E, ...
%                                         [], model.Mi.CS_z(E), ...
%                                         options_Opt_ROA);
                                    
% % [Origianl: CSQC + Valley QC]
% disp('===[ CS + Vallye ]===')
% Mi_all = func_stackingQCs(model.Mi,'all','ieQC',E);
% 
% [R0max_CSVal, info_CSVal, R0s_CSVal] = SDPOpt_ROAEstimation(model, E, ...
%  										[], Mi_all, ...
%                                         options_Opt_ROA);

%% Shape iteration
CS_eQC = {};
CS_ieQC = {'CS_z'};

CSValley_eQC = {};
CSValley_ieQC = {'CS_z';'Valley_z'};


% setting
niter = 3;
options_Shape.Opt_ROA = options_Opt_ROA;
options_Shape.verbose = true;

% disp('===[ CS ShapeIter ]===');
% [r_CS, info_CS, rs_CS, infos_CS] = func_ShapeIteration(model, ...
%                                                         CS_eQC, CS_ieQC, ...
%                                                         E0, niter, ...
%                                                         options_Shape);
%                        
% disp('===[ CS NewAlg ]===');                                                    
% [r_new, info_new] = func_ROAnewAlgorithm(model, CS_eQC, CS_ieQC, E0, options_Shape);                                                
        
% 
% disp('===[ CS_valley ShapeIter ]===');
% [r_CSValley, info_CSValley, rs_CSValley, infos_CSValley] = func_ShapeIteration(model, ...
%                                                         CSValley_eQC, CSValley_ieQC, ...
%                                                         E0, niter, ...
%                                                         options_Shape);                                             

disp('===[ CS_valley NewAlg ]===');  
options_Shape.alp_bisection = false;
[r_valley_new, info_valley_new] = func_ROAnewAlgorithm(model, CSValley_eQC, CSValley_ieQC, E0, options_Shape);                                                
            
%% save data

% save('ROA_2State');
save(['ROA_2State\', datestr(now, 'yyyymmdd_HHMM')]);

%% plotting

plotting_2State;
