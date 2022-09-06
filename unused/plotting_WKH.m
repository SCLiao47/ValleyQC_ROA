
data = zeros(6,3);

data(1,:) = [0.221, 0.0480, 0.0215];        % DAL
data(2,:) = [0.0383, 0.0063, 0.0021];       % Liu & Gayme
data(3,:) = [0.120, 0.0282, 0.0121];        % toso


data(4,:) = [0.1013, 0.0228, 0.0097];       % Kalur L-CSS

% data(5,:) = [0.1019, 0.0254, 0.0118];       % our 20211215    (loss
% data(6,:) = [0.1054, 0.0266, 0.0119];       % our 20220314

data(5,:) = r_data(1,:);
data(6,:) = r_data(2,:);                    % this setting


f1 = figure(1); clf

hand_lines = cell(6,1);
color_list = lines(6);

for i = 1:6
    hand_lines{i} = semilogy([5,10,15], data(i,:),'*-','color',color_list(i,:));    hold on;
end

% for i = 4:5
%     hand_lines{i} = semilogy(grid_Re, data(i,1),'*-','color',color_list(i,:)); hold on;
% end


set(findobj(gcf,'type','line'),'linewidth',2);
set(hand_lines{1},'color','k','linestyle','--');
% set(hand_lines{3},'linestyle','--');
% set(hand_lines{4},'linestyle','--');
% set(hand_lines{5},'linestyle','--');
set(hand_lines{6},'linewidth',3,'linestyle','-');   


grid on;

title('WKH Model B&T Parameter');
hl = legend([hand_lines{:}], {'Upper bound','Liu & Gayme','Toso','Kalur LCSS','Our work (20220314)', 'New Algorithm'});
set(hl,'fontsize',16);
set(hl,'location','southwest')

xlabel('\it Re','fontsize',16)
ylabel('$r^*$','fontsize',20,'interpreter','latex')
xlim([5,15]);
% ylim([0.001, 0.221])

% ylim([0.008, 0.221])

set(f1,'position',[865 292 745 433]);


% print('figure/WKH/compare', '-dpng', '-r600')