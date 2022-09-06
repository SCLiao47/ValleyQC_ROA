% Scripts that implements the largest spherical ROA estimation algorithm [1] for a 
% 2-state model in [2]. This analysis includes our proposed monomial formulation and valley QCs


% [Reference]
%   [1] A. Kalur, T. Mushtaq, P. Seiler, and M. S. Hemati, “Estimating regions
% 	of attraction for transitional flows using quadratic constraints,” IEEE
% 	Control Systems Letters, 2021.
%   [2] F. Amato, C. Cosentino, and A. Merola, “On the region of asymptotic
% 	stability of nonlinear quadratic systems,” in 2006 14th Mediterranean
% 	Conference on Control and Automation, pp. 1–5, IEEE, 2006.

clc; clear; close all;

%% models and settins
createModel = @model_WKH;
createModel_mon = @model_WKH_monomial;

% ===[ B&T ]===
which_model = 'WKH_B&T';
BnTParams = struct('ParamSet','B&T','lambda',1,'mu',1,'nu',1,'rho',1,'delta',1,'gamma',1);
modelParams = BnTParams;
% grid_Re = [5,10,15,20,25];  % grid of Reynold numbers to perform analysis
grid_Re = [5,10,15];
% grid_Re = [20,25];
% grid_Re = [5,10];
% grid_Re = [15];

NRe = length(grid_Re);

% test model
modelParams.Re = grid_Re(1);
model = createModel(modelParams);
model_mon = createModel_mon(modelParams);


%% Largest spherical ROA estimation
options_Opt_ROA = func_getOptions_SDP_ROA(1,-3,200,-6,false);

% options_Opt_ROA.Nalp = 20;
% options_Opt_ROA.verbose = true;

% create model
% model = createModel();

%% Checking implementation to match Liu & Gayme
% 
% % using spherical shape
% E0 = eye(model.nx);
% 
% % [Original: CSQC]
% disp('===[ CS ]===');
% [R0max_CS, info_CS, R0s_CS] = SDPOpt_ROAEstimation(model, E,...
%                                         [], model.Mi.CS_z(E), ...
%                                         options_Opt_ROA);
% 
% % [Origianl: CSQC + Valley QC]
% disp('===[ CS + Vallye ]===')
% Mi_all = func_stackingQCs(model.Mi,'all','ieQC',E);
% 
% [R0max_CSVal, info_CSVal, R0s_CSVal] = SDPOpt_ROAEstimation(model, E, ...
%  										[], Mi_all, ...
%                                         options_Opt_ROA);

%% Shape iteration
CS_eQC = {'lossless'};
CS_ieQC = {'CS_z'};

% Mon_eQC = {'lossless'};
Mon_eQC = {};

% Mon_ieQC = {'CS_z';'CS_phi';'Valley2_z'};         %20211215
% Mon_ieQC = {'CS_z';'CS_phi';'Valley2_z';'Valley2_phi';'Valley3_phi';'CrossP_z';'CrossP_neg'};
% Mon_ieQC = {'CS_z';'CS_phi';'Valley2_z';'Valley2_phi';'Valley3_phi';'CrossP_z'};    
Mon_ieQC = {'CS_z';'CS_phi';'Valley2_z';'Valley2_phi'};


% shaper iteration setting
niter = 3;
options_Shape.Opt_ROA = options_Opt_ROA;
options_Shape.verbose = true;
% options_Shape.verbose = false;

% initialize storage 
r_data = zeros(2,NRe);
info_data = cell(2,NRe);
rs_data = cell(2,NRe);
infos_data = cell(2,NRe);

for idx_Re = 1:NRe
    disp(['==== [ Re = ',num2str(grid_Re(idx_Re)),' ] ====']);
    
    modelParams.Re = grid_Re(idx_Re);
    
    % create model
    model = createModel(modelParams);
    model_mon = createModel_mon(modelParams);
    
    % initialize shape
    E0 = eye(model.nx);

%     disp('===[ L-CSS ]===');
%     [r_CS, info_CS, rs_CS, infos_CS] = func_ShapeIteration(model, ...
%                                                             CS_eQC, CS_ieQC, ...
%                                                             E0, niter,...
%                                                             options_Shape);
    disp('===[ Mon ]===');
    [r_Shape, info_Shape, rs_Shape, infos_Shape] = func_ShapeIteration(model_mon, ...
                                                            Mon_eQC, Mon_ieQC, ...
                                                            E0, niter, ...
                                                            options_Shape);
    r_data(1,idx_Re) = r_Shape;
    info_data{1,idx_Re} = info_Shape;
    rs_data{1,idx_Re} = rs_Shape;
    infos_data{1,idx_Re} = infos_Shape;
                                                        

  
    disp('===[ New Algorithm ]===');
    options_Shape.alp_bisection = true;
    [r_mon, info_mon] = func_ROAnewAlgorithm(model_mon, Mon_eQC, Mon_ieQC, ...
                                                E0, options_Shape);                                           
    rs_mon = [];
    infos_mon = [];                                            
                                                        
    r_data(2,idx_Re) = r_mon;
    info_data{2,idx_Re} = info_mon;
    rs_data{2,idx_Re} = rs_mon;
    infos_data{2,idx_Re} = infos_mon;
    
%     disp('===[ Alpha Gridding ]===');
%     options_Shape.alp_bisection = false;
%     [r_grid, info_grid] = func_ROAnewAlgorithm(model_mon, Mon_eQC, Mon_ieQC, ...
%                                                 E0, options_Shape);
%     rs_grid = [];
%     infos_grid = []; 
%     
%     r_data(1,idx_Re) = r_grid;
%     info_data{1,idx_Re} = info_grid;
%     rs_data{1,idx_Re} = rs_grid;
%     infos_data{1,idx_Re} = infos_grid;
end

%% hard-code data

plotting_WKH;