function make_summary_movie

load('img_dstruct.mat')

whitebg('white');
close all

try
    delete('summary_movie.mp4')
catch
    disp('no previous movie found')
end

mkdir('movie_frames')
cd('movie_frames')


[xymin, xymax] = determine_caxis(datastruct.pcd.raw_mstack, .4, .01);
[zymin, zymax] = determine_caxis(datastruct.pcd.raw_zstack, .65, .07);


pframe = 0;

zframe_num = 1;

wh = waitbar(0,'generating frames...');

hwc=findobj(wh,'Type','Patch');
set(hwc,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]) 

for frame_num = 1:10:size(datastruct.pcd.raw_mstack,3)-10;
%for frame_num = 1

close all
f1 = figure('color', 'w', 'units', 'normalized',...
    'position', [.01 .05 .95 .95], 'visible', 'off');

s1 = subplot(3, 5, [1:2 6:7]);
imagesc(datastruct.pcd.raw_mstack(:,:,frame_num));
set(gca, 'YDir', 'normal');
caxis([xymin xymax])
colormap(kjetsmooth(2^16));
axis equal tight off
title('x-y projection', 'fontsize', 25)

s2 = subplot(3, 5, [11 12]);


imagesc(rot90(mean(datastruct.pcd.raw_zstack(:,:,frame_num:frame_num+10), 3)));
caxis([zymin zymax])
colormap(kjetsmooth(2^16));
set(gca, 'YDir', 'normal');
daspect([60 10 1])
axis off tight
title('x-z projection', 'fontsize', 25)


s3 = subplot(3, 5, [3:5]);
plot(datastruct.pcd.tstamp_mstack, datastruct.pcd.roi_trace_dF(1,:))
hold on
plot([-1000 10000], [0 0], 'k')
z1 = scatter(datastruct.pcd.tstamp_mstack(frame_num), ...
    datastruct.pcd.roi_trace_dF(1, frame_num), 100);
set(z1, 'markeredgecolor', 'r', 'markerfacecolor', 'r')

xlim([0 datastruct.pcd.tstamp_mstack(end)]);
ylabel('dF/F', 'fontsize', 25)
box off


s4 = subplot(3,5, [8:10]);
plot(datastruct.pcd.tstamp_mstack, datastruct.pcd.set_temps, 'b')
hold on
xlim([0 datastruct.pcd.tstamp_mstack(end)]);
z1 = scatter(datastruct.pcd.tstamp_mstack(frame_num), ...
    datastruct.pcd.set_temps(frame_num), 100);
set(z1, 'markeredgecolor', 'r', 'markerfacecolor', 'r')
ylabel('set temps (°C)', 'fontsize', 25)
ylim([22 36])
box off

s5 = subplot(3,5, [13:15]);
plot(datastruct.pcd.tstamp_mstack, datastruct.pcd.measured_temps, 'b')
hold on
xlim([0 datastruct.pcd.tstamp_mstack(end)])
xlabel('time (sec)', 'fontsize', 25)

z1 = scatter(datastruct.pcd.tstamp_mstack(frame_num), ...
    datastruct.pcd.measured_temps(frame_num), 100);
set(z1, 'markeredgecolor', 'r', 'markerfacecolor', 'r')
ylabel('measured temps (°C)', 'fontsize', 25)

ylim([22 36])
box off

clear z1

pframe = pframe+1;

export_fig(['frame_' num2str(pframe, '%06d') '.bmp'], '-bmp', '-nocrop')


waitbar(frame_num / size(datastruct.pcd.raw_mstack,3))


end

close(wh)

disp('creating movie with ffmpeg')

success = system('ffmpeg -framerate 10 -i frame_%06d.bmp -c:v libx264 -r 30 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" ../summary_movie.mp4 ');
if success ~=0 
    path1 = getenv('PATH');
    path1 = [path1 ':/usr/local/bin'];
    setenv('PATH', path1);
    success = system('ffmpeg -framerate 10 -i frame_%06d.bmp -c:v libx264 -r 30 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" ../summary_movie.mp4 ');
end

cd('..')
rmdir('movie_frames', 's')

end

