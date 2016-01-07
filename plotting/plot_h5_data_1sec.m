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

f1 = figure('Color', 'w', 'units', 'normalized', 'position', [.01 .01 .9 .7]);
%whitebg('white')

hold on
xvals_ao = linspace(0,10,length(clk_data.recorded_data(:,1)));
xvals_ai = linspace(0,10,length(LED_signal));

LED_sec = clk_data.recorded_data(:,2)';
LED_sec = LED_sec(1:10000);
plot(xvals_ai(1:30000), diff(frame_counter(1:30001))+1.5, 'r', 'linewidth', 1.5)
plot(xvals_ai(1:30000), [diff(frame_out(1:30001))-1.5], 'k', 'linewidth', 1.5)

xlim([-.1 1.1])
ylim([-.1 5])

set(gca, 'XTick', [0 .5 1], 'YTick', [],...
       'FontSize', 22)
  

title('frame acquisition timestamps', 'FontSize', 25)
text(0, 4.5, 'frame out (ThorSync)', 'color', 'r', 'FontSize', 15)
text(0, 4.2, 'frame counter (ThorSync)', 'color', 'k', 'FontSize', 15)
text(0, 3.9, 'frame timestamps (ThorImage)', 'color', 'b', 'FontSize', 15)

img_samples = img_struct.tstamps(img_struct.tstamps<1);
z1 = scatter(img_samples, ones(1,length(img_samples))+.3, 200, 's');
set(z1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')

%ylabel(['voltage' char(0) '(relative)'])
xlabel('time(sec)')

cd(save_dir)
cd('plots')
export_fig('frame_stamp.pdf', '-pdf')
