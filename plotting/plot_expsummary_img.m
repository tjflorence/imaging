
h5_file = 'sync_signals.h5';
xml_file = 'file_locations.xml';

xy_project_b = 'ministacks_medfilt+ds_max.b';
xz_project_b = 'zxmip_medfilt+ds_max.b';

frame_xy_x = h5read(h5_file, '/settings/zmip/frame_dim_x');
frame_xy_y = h5read(h5_file, '/settings/zmip/frame_dim_y');

frame_xz_x = h5read(h5_file, '/settings/xmip/frame_dim_z');
frame_xz_y = h5read(h5_file, '/settings/xmip/frame_dim_x');

%% read out, reshape xy mip
f_xy = fopen(xy_project_b);

raw_out = fread(f_xy, 'uint16');
num_frames = numel(raw_out)/(frame_xy_x*frame_xy_y);
xy_project_data = reshape(raw_out, [frame_xy_y frame_xy_x num_frames]);

fclose(f_xy);

%% read out, reshape xz mip
f_xz = fopen(xz_project_b);

raw_out = fread(f_xz, 'uint16');
xz_project_data = reshape(raw_out, [frame_xz_x frame_xz_y num_frames]);

fclose(f_xz);


whitebg('black');
close all;

f1 = figure('color', 'k');

s1 = subplot(4,3,1:9)
imagesc(mean(xy_project_data(3:end-2, 3:end-2, :),3))
colormap(kjetsmooth(2^16))
title('x-y projection', 'fontsize', 25)
set(gca, 'YDir', 'normal')
axis equal off tight

s2 = subplot(4,3,11);
imagesc(rot90(mean(xz_project_data,3)))
colormap(kjetsmooth(2^16))
set(gca, 'YDir', 'normal');
title('x-z projection', 'fontsize', 25)
axis equal off tight
daspect([10 2 1])

p2 = get(s2, 'Position')
set(s2, 'Position', [p2(1)-.75 p2(2)-.01 p2(3)+1.5 p2(4)])

p1 = get(s1, 'Position')
set(s1, 'Position', p1*.9)

p1a = get(s1, 'Position')
set(s1, 'Position', [p1a(1)+.05 p1a(2)+.05 p1a(3) p1a(4)])
export_fig('summary_projection.pdf', '-pdf', '-zbuffer')