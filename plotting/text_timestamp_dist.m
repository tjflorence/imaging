close all

save_dir = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/img_acq_001';

h5_datafile = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/sync001/Episode001.h5';
clk_file    = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/clk_period_10_sec_10_hz_001.mat';
img_file    = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/img_acq_001/img_struct.mat';

LED_signal  = h5read(h5_datafile, '/AI/LED');
clk_signal  = h5read(h5_datafile, '/AI/clk');
trigger     = h5read(h5_datafile, '/AI/ext trigger');

frame_counter   = h5read(h5_datafile, '/CI/Frame Counter');
frame_in        = h5read(h5_datafile, '/DI/Frame In');
frame_out       = h5read(h5_datafile, '/DI/Frame Out');


load(img_file);
load(clk_file);

 cmap = [204, 54, 19;...
        229, 39, 239;...
        31, 255, 253;...
        255, 240, 18;...
        9 243 27]/255;

%% hist of pixel progression
f1 = figure('Color', 'k', 'units', 'normalized',...
    'Position', [0.1300    0.1100    0.7750    0.8150])

hist(diff(img_struct.tstamps),25);
set(get(gca,'child'),'FaceColor', cmap(3,:),'EdgeColor','k');
hold on

%plot([0 0], [0 100], 'color', 'w', 'linewidth', 2)
%plot([0 5], [0 0], 'color', 'w', 'linewidth', 2)

box off


%text(3, 70, ['expected frame' char(10) 'progression'],...
 %   'fontsize', 15, 'fontweight', 'bold')
xlim([0 .16])
ylim([0 20])

set(gca, 'XTick', [.05 .1 .15], 'YTick', [5 10 15], 'FontSize', 22)

xlabel(['inter-frame interval' char(10) '(sec)' ], 'FontSize', 25)
ylabel('count', 'FontSize', 25)

cd(save_dir)
export_fig('tstamp_distribution.pdf', '-pdf')
