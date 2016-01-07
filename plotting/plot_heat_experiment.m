homedir = pwd;


heat_files = dir('env_heatmap*');
load('roi_data.mat');
load('img_dstruct.mat')

filenum = 1;
roi_num = 1;

load(heat_files(1).name);

for roi_num = 1:length(roi_struct)
whitebg('black');
close all;

f1 = figure('color', 'k', 'units', 'normalized', 'Position', [.028 0 .4923 .9048]);

s1 = subplot(5,5, [1:2 6:7]);
imagesc(mean(datastruct.pcd.corr_mstack, 3));
set(gca, 'YDir', 'normal')
axis equal off tight
caxis([60 80])
freezeColors();
hold on
title('summary image', 'fontsize', 30, 'Fontweight', 'bold');


s2 = subplot(5,5, [4:5 9:10]);
whitebg('black')
imagesc(mean(datastruct.pcd.corr_mstack, 3));
set(gca, 'YDir', 'normal')
colormap([0,0,0; 0,0,0; gray(8)])
caxis([60 80])
axis equal off tight
hold on
title(['ROI ' num2str(roi_num)], 'fontsize', 30, 'fontweight', 'bold',...
    'color', roi_struct(roi_num).cmap);


scat_h = fill(roi_struct(roi_num).xy(:,1), roi_struct(roi_num).xy(:,2),  roi_struct(roi_num).cmap);
  set(scat_h, 'MarkerFaceColor', roi_struct(roi_num).cmap, 'MarkerEdgeColor', 'w');
  

%% plot trial 1
s3 = subplot(5,5, [11:15]);
start_x = expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame(1));
end_x = expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame(end));

plot([-1000 1000], [0 0], 'color', 'w');
hold on

z1 = fill([240 480 480 240], [-.5 -.5 .5 .5], [.4 .4 .4]);
set(z1, 'facealpha', .5)

z1 = fill([480 720 720 480], [-.5 -.5 .5 .5], [.7 .7 .7]);
set(z1, 'facealpha', .5)

for aa = 1:6
    z1 = fill([60 120 120 60]+((aa-1)*120), [-.5 -.5 .5 .5], [78 166 255]./255);
    set(z1, 'facealpha', .5)
end

plot(expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame),...
    expr.c_trial.data.img.roi_trace_dF(roi_num,:));

ymin = min(expr.c_trial.data.img.roi_trace_dF(roi_num,:))-.1;
ymax = max(expr.c_trial.data.img.roi_trace_dF(roi_num,:))+.1;
xlim([start_x end_x])
ylim([ymin ymax])
box off

text(0, ymax+.05, 'dark', 'fontsize', 20)
text(240, ymax+.05, 'light', 'fontsize', 20)
text(480, ymax+.05, 'VR img', 'fontsize', 20)

set(gca, 'XTick', [0:120:720], 'XTickLabel', {}, 'FontSize', 20)

%% plot trial 2
load(heat_files(2).name)

s4 = subplot(5,5, [16:20]);
start_x = expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame(1));
end_x = expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame(end));

plot([-1000 1000], [0 0], 'color', 'w');
hold on

z1 = fill([240 480 480 240], [-.5 -.5 .5 .5], [.4 .4 .4]);
set(z1, 'facealpha', .5)

z1 = fill([480 720 720 480], [-.5 -.5 .5 .5], [.7 .7 .7]);
set(z1, 'facealpha', .5)

for aa = 1:6
    z1 = fill([60 120 120 60]+((aa-1)*120), [-.5 -.5 .5 .5], [78 166 255]./255);
    set(z1, 'facealpha', .5)
end

plot(expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame),...
    expr.c_trial.data.img.roi_trace_dF(roi_num,:));

xlim([start_x end_x])
ylim([ymin ymax])
box off

set(gca, 'XTick', [0:120:720], 'XTickLabel', {} ,'FontSize', 20)
ylabel('dF/F')

%% plot trial 3
load(heat_files(3).name)

s4 = subplot(5,5, [21:25]);
start_x = expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame(1));
end_x = expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame(end));

plot([-1000 1000], [0 0], 'color', 'w');
hold on

z1 = fill([240 480 480 240], [-.5 -.5 .5 .5], [.4 .4 .4]);
set(z1, 'facealpha', .5)

z1 = fill([480 720 720 480], [-.5 -.5 .5 .5], [.7 .7 .7]);
set(z1, 'facealpha', .5)

for aa = 1:6
    z1 = fill([60 120 120 60]+((aa-1)*120), [-.5 -.5 .5 .5], [78 166 255]./255);
    set(z1, 'facealpha', .5)
end

plot(expr.c_trial.data.timestamp(expr.c_trial.data.img.trial_frame),...
    expr.c_trial.data.img.roi_trace_dF(roi_num,:));

xlim([start_x end_x])
ylim([ymin ymax])
box off

set(gca, 'XTick', [0:120:720] ,'FontSize', 20)
xlabel('time (sec)')


cd('plots')
export_fig(['heatpulse_ROI_' num2str(roi_num, '%02d') '.pdf'], '-pdf', '-zbuffer')
cd(homedir)


end
