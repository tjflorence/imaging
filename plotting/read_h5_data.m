save_dir = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/img_acq_001'

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

f1 = figure('Color', 'k', 'units', 'normalized', 'position', [.01 .01 .9 .7]);
whitebg('black')

s1 = subplot(2,1,1);

set(s1, 'Color', 'w');
hold on
xvals = linspace(0,10,length(clk_data.recorded_data(:,1)));
plot(xvals, clk_data.recorded_data(:,1)'+6, 'k')
plot(xvals, clk_data.recorded_data(:,2)', 'r')
plot(xvals, clk_data.recorded_data(:,3)'-6, 'b')

xlim([-.1 10.1])
ylim([-15 15])

set(gca, 'XTick', [], 'YTick', [-15 -5 5 15],...
       'YTickLabel', {'0', '10', '20', '30'}, 'FontSize', 22)
   

title('command computer', 'FontSize', 25)
text(0, 13.5, 'trigger', 'color', 'k', 'FontSize', 15)
text(1, 13.5, 'LED', 'color', 'r', 'FontSize', 15)
text(2, 13.5, 'clock', 'color', 'b', 'FontSize', 15)

s2 = subplot(2,1,2);
set(s2, 'Color', 'k')

xvals = linspace(0,10, length(LED_signal));
hold on

plot(xvals, trigger+10, 'color', 'w')
plot(xvals, LED_signal, 'color', cmap(1,:))
plot(xvals, clk_signal-6, 'color', cmap(3,:))

set(gca, 'XTick', [2 4 6 8], 'YTick', [-15 -5 5 15],...
       'YTickLabel', {'0', '10', '20', '30'}, 'FontSize', 22)

xlim([-.1 10.1])
ylim([-15 15])

yl_h = ylabel(['Voltage' char(10) '(relative)']);
xlabel('time (sec)')
ylabel_pos = get(yl_h, 'Position');
set(yl_h, 'Position', [ylabel_pos(1) ylabel_pos(2)+20 ylabel_pos(3)]);
title('imaging computer', 'FontSize', 25)
text(0, 13.5, 'trigger', 'color', 'w', 'FontSize', 15, 'fontweight', 'bold')
text(1, 13.5, 'LED', 'color', cmap(1,:), 'FontSize', 15, 'fontweight', 'bold')
text(2, 13.5, 'clock', 'color', cmap(3,:), 'FontSize', 15, 'fontweight', 'bold')

cd(save_dir)
cd('plots')
export_fig('signals_10_sec.pdf', '-pdf')
