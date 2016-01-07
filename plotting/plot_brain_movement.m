
function plot_brain_movement()

homedir = pwd;

load('img_dstruct.mat')
load('roi_data.mat')

for jj = 1:length(roi_struct)
    
xvals = 2*roi_struct(jj).xy(:,1);

select_ind = round(mean(xvals));

zmat = nan(size(datastruct.pcd.raw_zstack,2), ceil(size(datastruct.pcd.raw_zstack,3)/10) );
ind_vec = nan(1, ceil(size(datastruct.pcd.raw_zstack,3)/10)-1);

p = 0;

for ii = 1:10:(size(datastruct.pcd.raw_zstack,3)-10);
    
    p = p+1;
    
    c_frame = rot90(mean(datastruct.pcd.raw_zstack(:,:,ii:ii+10), 3));
    c_col = c_frame(:,select_ind);
    
    zmat(:,p) = c_col;
    ind_vec(p) = ii;
    
end



tstamps = datastruct.pcd.tstamp_mstack(ind_vec);
mtemps = datastruct.pcd.measured_temps(ind_vec);
stemps = datastruct.pcd.set_temps(ind_vec);
dF = datastruct.pcd.roi_trace_dF(jj,:);
dF = dF(ind_vec);

f1 = figure('Position', [23 1052 1237 1502], 'visible', 'off');

s1 = subplot(3,1,1);
colormap(kjetsmooth(2^16));
[cmin, cmax] = determine_caxis(datastruct.pcd.raw_zstack, .65, .07);
imagesc(zmat)
caxis([cmin cmax])
set(gca,  'YDir', 'normal')
box off
axis off

s2 = subplot(3,1,2);
plot([tstamps(1) tstamps(end)], [0 0], 'k')
hold on
plot(tstamps, dF, 'linewidth', 1.5, 'color', roi_struct(jj).cmap)

xlim([tstamps(1) tstamps(end)])
box off
set(gca, 'XTickLabel', {}, 'FontSize', 25)
ylabel('dF/F', 'fontsize', 25)

s3 = subplot(3,1,3);
hold on
plot(tstamps, stemps, 'k', 'linewidth', 1)
plot(tstamps, mtemps, 'r', 'linewidth', 1.2);
set(gca, 'FontSize', 25)
ylim([22 37])

xlim([tstamps(1) tstamps(end)])
box off

xlabel('time (sec)', 'fontsize', 30)
ylabel('temp (°C)')

set(gcf, 'units', 'inches')
pos = get(gcf, 'Position');
set(gcf, 'paperpositionmode', 'manual', 'paperposition', [0 0 pos(3)*.5 pos(4)*.5]);

cd('plots')
print(['zmovement_summary_ROI_' num2str(jj, '%02d') '.pdf'], '-dpdf')

cd(homedir)

close all
end


end