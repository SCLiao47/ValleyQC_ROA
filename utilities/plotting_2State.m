%% plot phase portrait
model = model_2State();

Tspan = 0.1;
Np = 10;

figure(1); clf;

h_max = 6;
h_min = -4;
v_max = 4;
v_min = -4;

% trajectories from boundary
Tspan = 0.4;
Np = 10;

Vgrid = linspace(v_min,v_max, Np)';
Hgrid = linspace(h_min,h_max, Np)';
Inp = ones(Np,1);

sideD = [Hgrid v_min*Inp];      % down side
sideR = [h_max*Inp Vgrid];      % right side
sideU = [Hgrid v_max*Inp];      % upper side
sideL = [h_min*Inp Vgrid];      % left side

% Trimming
sideD([2,3,6,7,8,9],:) = [];
sideR(7,:) = [];
sideU([1,3,5,6,7:end],:) = [];
sideL([2,4,6,8,9,10],:) = [];

% adding
sideR = [sideR;
         6, 0.9;
         6, 1.7];

x0_grid = [sideD; sideR; sideU; sideL];

hp = [];
for i = 1:size(x0_grid,1)
    x0 = x0_grid(i,:);
    [T,X] = ode45(@(t,x) model.func_ode(t,x), [0, Tspan], x0);
    
    % plot traj
    hp = [hp; plot(X(:,1),X(:,2),'color',[1,1,1]*0.5)]; hold on;
end
set(hp,'linewidth',2);

% unstable manifold
[T,X] = ode45(@(t,x) model.func_ode(t,x), [0, Tspan], [9.8600473 0]);
plot(X(:,1),X(:,2),'r--','linewidth',2.5); 

% add vector field
[X_mesh,Y_mesh] = meshgrid(h_min:0.5:h_max, v_min:0.5:v_max);
U_mesh = 0*X_mesh;
V_mesh = 0*X_mesh;

for i = 1:size(X_mesh,1)
    for j = 1:size(X_mesh,2)
        x = [X_mesh(i,j); Y_mesh(i,j)];
        xdot = model.func_ode(1,x);
        
        U_mesh(i,j) = xdot(1);
        V_mesh(i,j) = xdot(2);
    end
end

hold on
quiver(X_mesh,Y_mesh,U_mesh,V_mesh,'color',[1,1,1]*0.5, 'LineWidth',1)
% hold off

grid on

axis equal 
xlim([-4 6])
ylim([-4 4])

xlabel('x_1','fontsize',16);
ylabel('x_2','fontsize',16);

% print('figure/MCC/PhasePortrait', '-dpng', '-r600')

%% plot ROA
Inx = eye(2);
cxy = [0,0];

% Set containment condition: 
%   x'x <= x'Px <= q*x'x = c 

% from data
% r_CS = 2.7355;
% r_CSValley = 3.5224;

% ROA sphere: x'*x = R0^2
% h_ROA_mon_CS = plot_ellipse(2.7349, cxy, Inx, [0 0.4470 0.7410]);
h_ROA_mon_CS = plot_ellipse(r_CS, cxy, Inx, [0 0.4470 0.7410]);
set(h_ROA_mon_CS,'linewidth',3);

% h_ROA_mon_CSValley = plot_ellipse(3.5706, cxy, Inx, [0.4660 0.6740 0.1880]);
h_ROA_mon_CSValley = plot_ellipse(r_CSValley, cxy, Inx, [0.4660 0.6740 0.1880]);
set(h_ROA_mon_CSValley,'linewidth',3);
hs = [h_ROA_mon_CS; h_ROA_mon_CSValley];

%% for paper
axis equal 
xlim([-4 6])
ylim([-4 4])

% set(gcf,'position',[242 253 691 549]);
hl = legend([h_ROA_mon_CS; h_ROA_mon_CSValley], {'Set 1', 'Set 2'});

set(hl,'fontsize',16,'location','southeast');

xlabel('x_1','fontsize',16);
ylabel('x_2','fontsize',16);

print('Data/figure/ROA_2State', '-dpng', '-r600')
