
%% set experiment location
experiment_dir = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/img_acq_001';


%% set experiment parameters
flicker_freq = 10;
actual_img_freq = 10.3;

img_struct = load_tseries(experiment_dir);


% color
 cmap = [204, 54, 19;...
        229, 39, 239;...
        31, 255, 253;...
        255, 240, 18;...
        9 243 27]/255;

%% set column range to use for analysis
n_range = [560 700];

data_vals = img_struct.img(:, n_range(1):n_range(2),:);

mean_column = mean(data_vals,2);
dims = size(mean_column);

mean_2d = reshape(mean_column, [dims(1) dims(3)]);

nwise_diff = diff(mean_2d);
blank_mat = zeros(size(nwise_diff));
blank_mat(nwise_diff<-400) = -1;
blank_mat(nwise_diff>300) = 1;

[yon xon] = find(blank_mat==1);
[yoff xoff] = find(blank_mat==-1);


yon_diff = diff(yon);

yon_to_plot = yon;
yon_data = yon_diff;
yon_to_plot([yon_diff<-700]) = nan;
yon_data(yon_diff<-700) = yon_diff(yon_diff<-700)+size(blank_mat,1);

whitebg('black')

close all

cd(experiment_dir)
mkdir('plots')
cd('plots')


f1 = figure('Color', 'k', 'units', 'normalized', 'Position', [.1 .1 .8 .5])

imagesc(mean_2d)
colormap(gray)
set(gca, 'Ydir', 'normal')
box off

xlabel('frame', 'fontsize' ,25)
ylabel(['LED saturation ' char(10) 'extent'], 'fontsize', 25)

set(gca, 'XTick', [30 60 90], 'YTick', [200 400 600], 'FontSize', 22)
export_fig('phase_progression_raw.png', '-png')


hold on
plot(xon, yon_to_plot, 'color', cmap(2,:), 'linewidth', 5)
drawnow
export_fig('phase_progression_on_pix_located.pdf', '-pdf')

hold off

close all

actual_img_dur = 1/actual_img_freq;
row_acq_dur = actual_img_dur/size(blank_mat,1);

expected_progression = abs((1/flicker_freq)-(1/actual_img_freq));

%% hist of pixel progression
f1 = figure('Color', 'k', 'units', 'normalized',...
    'Position', [0.1300    0.1100    0.7750    0.8150])

yon_data_sec = yon_data*row_acq_dur;
hist(yon_data_sec*1000,25);
set(get(gca,'child'),'FaceColor', cmap(3,:),'EdgeColor','k');
hold on

plot([0 0], [0 100], 'color', 'w', 'linewidth', 2)
plot([0 5], [0 0], 'color', 'w', 'linewidth', 2)

plot([expected_progression expected_progression]*1000, [-100 100], 'linewidth', 4)
box off

ylim([0 100])
xlabel(['observed frame progression' char(10) '(msec)' ], 'FontSize', 25)
ylabel('count', 'FontSize', 25)

text(3, 70, ['expected frame' char(10) 'progression'],...
    'fontsize', 15, 'fontweight', 'bold')
xlim([0 4])
ylim([0 100])

set(gca, 'XTick', [1 2 3], 'YTick', [30 60 90], 'FontSize', 22)

export_fig('progression_histogram.pdf', '-pdf')


%% make multicolor 2-gaussian distribution
f1 = figure('Color', 'k', 'units', 'normalized', 'Position', [0.1300    0.1100    0.7750    0.8150])

all_pix = mean_2d(:);
low_vals = all_pix(all_pix<600);
high_vals = all_pix(all_pix>600);

xrange = linspace(min(min(mean_2d)), max(max(mean_2d)), 500); 
[n1, xout1] = hist(low_vals,xrange);
bar(xout1,n1, 'g');
hold
[n2, xout2] = hist(high_vals,xrange);
bar(xout2,n2,'w');

box off

plot([5 5], [0 5000], 'Color', 'w', 'linewidth', 2)
plot([0 1500], [5 5], 'Color', 'w', 'linewidth', 2)

xlim([0 1200])
ylim([0 2700])

set(gca, 'XTick', [300 600 900], 'YTick', [], 'FontSize', 22)

xlabel('pixel value', 'fontsize', 25)
ylabel('count', 'fontsize', 25)

export_fig('pixel intensity distribution.pdf', '-pdf')

close all


%% convert progression to time, plot timestamp (thorsync) hist and 
% led-estimated hist together
whitebg('black')

% convert progression to actual frame timing
y_pix = size(img_struct.img, 1);
flashes_per_sec = flicker_freq;
yon_data(yon_diff<-700) = yon_diff(yon_diff<-700)+size(blank_mat,1);

pix_per_flash = yon_data+y_pix;
pix_per_sec = pix_per_flash*flashes_per_sec;
sec_per_pix = double(1./pix_per_sec);

sec_per_frame = y_pix*sec_per_pix;

yon_frame_time = (yon_data+y_pix)/pix_per_sec;

% get data from thorsync files
h5_datafile = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/sync001/Episode001.h5';
clk_file    = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/clk_period_10_sec_10_hz_001.mat';
img_file    = '/Users/florencet/Documents/matlab_root/imaging/meeting_and_analysis/201504/timing_tests/20150414/tsync_period_10_sec_10_hz/img_acq_001/img_struct.mat';

LED_signal  = h5read(h5_datafile, '/AI/LED');
clk_signal  = h5read(h5_datafile, '/AI/clk');
trigger     = h5read(h5_datafile, '/AI/ext trigger');

frame_counter   = h5read(h5_datafile, '/CI/Frame Counter');
frame_in        = h5read(h5_datafile, '/DI/Frame In');
frame_out       = h5read(h5_datafile, '/DI/Frame Out');

acq_hz = 30000;

frame_out   = diff(frame_out);
frame_inds  = find(frame_out>0);
inter_frame_frames = diff(frame_inds);

interframe_seconds = inter_frame_frames/acq_hz;

f1 = figure('Color', 'k', 'units', 'normalized',...
    'Position', [0.0649    0.3105    0.4857    0.5724])

xrange = linspace(min(min([sec_per_frame' interframe_seconds])),...
    max(max([sec_per_frame' interframe_seconds])), 50);


s1 = subplot(2,2,1:2)
[n1, xout1] = hist(interframe_seconds,xrange);
bar(xout1,n1, 'g');
xlim([0.0955 .1])
ylim([0 100])
ylabel('count', 'FontSize', 25)
title('interframe interval (ThorSync)', 'Color', 'g', 'FontSize', 25, 'FontWeight', 'bold')
set(gca, 'Ytick', [0 50 100], 'FontSize', 22)
box off

hold
s2 = subplot(2,2,3:4)
[n2, xout2] = hist(sec_per_frame,xrange);
bar(xout2,n2,'w');
xlim([0.0955 .1])
ylim([0 100])
ylabel('count', 'FontSize', 25)
xlabel('sec', 'FontSize', 25)
title('interframe interval (ThorImage/LED)', 'Color', 'w', 'FontSize', 25, 'FontWeight', 'bold')
set(gca, 'Ytick', [0 50 100], 'FontSize', 22)
box off

export_fig('thorsync_vs_thorimage.pdf', '-pdf')


hold off
%legend('First Series','Second Series') % add legend



close all

