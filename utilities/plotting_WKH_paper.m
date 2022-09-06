
load('ROA_WKH_BnT_20211215.mat');

data = zeros(5,3);

data(1,:) = [rs_data{1,1}(1), rs_data{1,2}(1), rs_data{1,3}(1)];    % Liu & Gayme
data(2,:) = r_data(1,:);                    % Kalur L-CSS
data(3,:) = [0.120, 0.0282, 0.0121];        % toso
data(4,:) = r_data(2,:);                    % our
data(5,:) = [0.221, 0.0480, 0.0215];        % DAL

f1 = figure(1); clf

lines = cell(5,1);



for i = 1:5
    lines{i} = semilogy(grid_Re, data(i,:),'*-');    hold on;
end

% lines{1} = semilogy(grid_Re, data(1,:),'*-','color',[0.8500 0.3250 0.0980]);       hold on;
% lines{2} = semilogy(grid_Re, data(2,:),'*-','color',[0 0.4470 0.7410]);         hold on;
% lines{3} = semilogy(grid_Re, data(3,:),'*-','color',[0.4940 0.1840 0.5560]);
% lines{4} = semilogy(grid_Re, data(4,:),'*-','color',[0.9290 0.6940 0.1250]);
% lines{5} = semilogy(grid_Re, data(5,:),'*-','color','k');
% 
% for i = 1:4
%     lines{i} = [];
% end
% lines{1} = [];
% lines{2} = [];
% lines{3} = [];
% lines{4} = [];
% lines{5} = [];


set(findobj(gcf,'type','line'),'linewidth',2);
set(lines{4},'linewidth',3,'linestyle','-');
set(lines{5},'linestyle','--');
grid on;

% title('Default');
% hl = legend([lines{:}], {'[1]','[2]','Our work','Upper bound'});
hl = legend([lines{:}], {' Direct sim.'});
% hl = legend([lines{:}], {' Kalur et al. [3]',' Our work',' Direct sim.'});
% hl = legend([lines{:}], {' Kalur et al. [3]',' Toso et al. [6]',' Our work',' Direct sim.'});
set(hl,'fontsize',16);

xlabel('\it Re','fontsize',16)
ylabel('$r^*$','fontsize',20,'interpreter','latex')
xlim([5,15]);
ylim([0.001, 0.221])

ylim([0.008, 0.221])

set(f1,'position',[865 292 745 433]);


print('figure/WKH/compare', '-dpng', '-r600')