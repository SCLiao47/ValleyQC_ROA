


model = model_CDC3States();

%% setting
r0 = 2.428391;     % initial radius
np = 1e4;   % number of point

% ode setting
odeopt = odeset('Events', @EventFcn);
Tspan = [0,20];

% looping setting
s = 0.001;
del = 0.000001;
tol = 1e-5;

flagCount = 1000;

%% unstable x0

% from eigenvector
[V,D] = eig(model.A);

rbar_eig = 4.818;
vbar_eig = V(:,2);

% from sampling
rbar_sam(1) = [2.4302];
vbar_sam(:,1) = [-0.9076; 0.3476; 0.2353];

rbar_sam(2) = [2.42834829];
vbar_sam(:,2) = [-0.9082; 0.3322; 0.2546];

%% simulation
flag = 0;

r = r0;

v_Diverge = [vbar_eig, vbar_sam];
r_Diverge = [rbar_eig, rbar_sam];
v_EdgeCase = [];
r_EdgeCase = [];

rbar = min(r_Diverge);

while flag < flagCount
    fprintf('\nr = %.6f, rbar = %.6f, flag = %i', r, rbar, flag);
    
    % reset flag
    isDiverge = false;
    
    %% sample initial condition on unit sphere with radius r
    % on cube [-1,1]
    x0s = 2*rand(3, np)-1;   
    
    % reject sample outside sphere
    x0sNorm = vecnorm(x0s);
    idx = or(x0sNorm > 1, x0sNorm < 0.001);
    x0s(:,idx) = [];
    
    % reject sample that is too samll
    
    % normalize the vectors
    x0s = x0s./vecnorm(x0s);
    
    % add an unstable direction found by eigenvector
    x0s = [x0s, v_Diverge];
    
    % scale the sphere to r
    x0s = r*x0s;
    
    %% all initial conditions
    for i = 1:size(x0s,2)
        x0 = x0s(:,i);
        
        % run ode45
        [T,X,te,xe,ie] = ode45(@(t,x) model.dyn(t,x), Tspan, x0, odeopt);
        
        if isempty(ie)
%             warning('something interesting happended!');
            ie = 3;
        end
            
        % check event flag
        switch ie
            case 1     % converge to 0 
                % do nothing
                
            case 2      % diverge
                isDiverge = true;
                
                v_Diverge = [v_Diverge, x0/r];
                r_Diverge = [r_Diverge, r];
                break;
                
            otherwise
                warning('Unforeseen even shows! Store to data and continue');
                v_EdgeCase = [v_EdgeCase, x0/r];
                r_EdgeCase = [r_EdgeCase, r];
        end
    end
    
    if isDiverge
        fprintf('Diverge!');
        
        flag = 0;
        
        rbar = r;
        
        r = r*(1 - s);
    else        
        if r*(1+s) > rbar
            r = min(rbar-del, r + del);
        else
            r = r*(1+s);
        end
        
        if abs(rbar - r)/rbar < tol
            flag = flag + 1;
        end
    end
end

%%
disp([]);
disp(['The simulation end. rbar = ', num2str(rbar)]);
v_Diverge(:,end)

%%
function [value,isterminal,direction] = EventFcn(t,x)
    % Event1: converging to origin
    value(1) = norm(x) - 0.7051;    % ROA estimated via CSQC
    isterminal(1) = 1;
    direction(1) = 0;
    
    % Event2: diverging
    value(2) = norm(x) - 100;
    isterminal(2) = 1;
    direction(2) = 0;
end
