function bin2struct(process_dir)

cd(process_dir);

xml_struct = xml2struct('file_locations.xml');
idir = xml_struct.file_dirs.idir.Text;
sync_dir = xml_struct.file_dirs.syncdir.Text;
save_dir = xml_struct.file_dirs.savedir.Text;
homedir = save_dir;
b_dir = homedir;

acq_hz = 10000;

%% read data from thorsync file
cd(sync_dir);
exp_h5_f = 'Episode001.h5';

pz_pos  = h5read(exp_h5_f, '/AI/Piezo Monitor');
clk     = h5read(exp_h5_f, '/AI/clk');
trial   = h5read(exp_h5_f, '/AI/trial');
frame_out = h5read(exp_h5_f, '/DI/Frame Out');
measured_temp = (h5read(exp_h5_f, '/AI/measured_temp')*3/2)+30; % these are flipped, need to fix
set_temp = (h5read(exp_h5_f, '/AI/set_temp')*3/2)+30;
sync_tstamps = (0:1:length(pz_pos))/acq_hz;

%% read pre-processed data
% load ministack binary file
cd(b_dir)
mini_f = 'ministacks_medfilt+ds_max.b';
z_f = 'zxmip_medfilt+ds_max.b';
stack_h5_f = 'sync_signals.h5';

% load from hdf5 file
laser_mW = h5read(stack_h5_f, '/settings/laser_power');
frame_dim_x = h5read(stack_h5_f, '/settings/zmip/frame_dim_x');
frame_dim_y = h5read(stack_h5_f, '/settings/zmip/frame_dim_y');

frame_xz_x = h5read(stack_h5_f, '/settings/xmip/frame_dim_z');
frame_xz_y = h5read(stack_h5_f, '/settings/xmip/frame_dim_x');

idx_vec = h5read(stack_h5_f, '/sync_data/idx_vec');
clk_mstack = h5read(stack_h5_f, '/sync_data/clk_vec');
mtemp_vec = h5read(stack_h5_f, '/sync_data/measured_temp');
stemp_vec = h5read(stack_h5_f, '/sync_data/set_temp');
trial_vec = h5read(stack_h5_f, '/sync_data/trial_vec');
framestack_id = h5read(stack_h5_f, '/sync_data/framestack_id');

disp('loading ministack file')
f_h = fopen(mini_f, 'r');
raw_out = fread(f_h, 'uint16');
num_frames = length(raw_out)/(frame_dim_x*frame_dim_y);
mstack_mat = reshape(raw_out, [frame_dim_y frame_dim_x num_frames]);
clear raw_out;
fclose(f_h);
disp('file loaded')

disp('loading zx file');
z_h = fopen(z_f, 'r');
raw_out = fread(z_h, 'uint16');
zstack_mat = reshape(raw_out, [frame_xz_x frame_xz_y num_frames]);
clear raw_out;
fclose(z_h);
disp('file loaded')

%% collect everything into a neat structure
% settings 
datastruct.settings.xy_frame_dim = [frame_dim_x frame_dim_y];
datastruct.settings.zx_frame_dim = [frame_xz_x frame_xz_y];
datastruct.settings.laser_mW = laser_mW;
datastruct.settings.num_mstacks = num_frames;

% ministack (processed) data
datastruct.pcd.raw_mstack = mstack_mat;
datastruct.pcd.raw_zstack = zstack_mat;
datastruct.pcd.tstamp_mstack = sync_tstamps(idx_vec);
datastruct.pcd.measured_temps = mtemp_vec;
datastruct.pcd.set_temps = stemp_vec;
datastruct.pcd.trial_vec = trial_vec;
datastruct.pcd.clk_mstack = clk_mstack;
datastruct.pcd.idx_vec = idx_vec;

% raw data from experiment
datastruct.raw.pz_pos   = pz_pos;
datastruct.raw.clk      = clk;
datastruct.raw.trial    = trial;
datastruct.raw.frame_out = frame_out;
datastruct.raw.timer_stamp  = sync_tstamps;
datastruct.raw.measured_temps = measured_temp;
datastruct.raw.set_temps = set_temp;
datastruct.raw.framestack_id = framestack_id;

%% now, create trial ids and count behavior simulation frames
datastruct.raw.trial_id = nan(1, length(datastruct.raw.trial));
datastruct.raw.trial_frame = nan(1, length(datastruct.raw.trial));

t_id = diff(datastruct.raw.trial);
t_starts = find(t_id>1);
t_ends = find(t_id<-1);

if length(t_starts) ~= length(t_ends)
    disp('number of trial starts and ends unequal, attempting to fix')
    t_start_fix = t_starts(1);
    t_diff = diff(t_starts);
    
    for zz = 1:length(t_diff)
        if t_diff(zz) > 10
            t_start_fix = [t_start_fix, t_starts(zz+1)];
        end
    end
        
    if length(t_start_fix)==length(t_ends)
           
        t_starts = t_start_fix;
        disp('fixed, motherfucka!')
        
    else
        
        disp('did not fix, but continuing')
    end
    
    
else
    disp('num trial starts == num trial ends')
end

disp('finding trial starts')
for ii = 1:length(t_ends)
   
    t_start_i = t_starts(ii);
    t_end_i = t_ends(ii);
    
    datastruct.raw.trial_id(t_start_i:t_end_i) = ii;
    
    
end
disp('trial starts found')


disp('counting behavioral frames')
wait_h = waitbar(0,'counting behavioral frames by trial');
hw=findobj(wait_h,'Type','Patch');
set(hw,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);

num_trials = max(datastruct.raw.trial_id);
for ii = 1:num_trials
    
    t_frames = find(datastruct.raw.trial_id==ii);
    clk_frames = datastruct.raw.clk(t_frames);
    d_clk = diff(clk_frames);
    d_clk(d_clk<-2) = .1;

    bframe_counter = zeros(1, length(d_clk));
    bframe_num = 1;

    for aa = 1:length(bframe_counter);
        tickval = d_clk(aa);
        if tickval > 0.07
            bframe_num = bframe_num+1;
        end
    
        bframe_counter(aa) = bframe_num;
    end

    datastruct.raw.trial_frame(t_frames(1:end-1)) = bframe_counter;
    wait_h = waitbar(   (ii) / num_trials);


end
disp('behavioral frames counted')
delete(wait_h)

datastruct.pcd.trial_id = nan(1, length(datastruct.pcd.tstamp_mstack));
datastruct.pcd.frame_id = nan(1, length(datastruct.pcd.tstamp_mstack));
disp('matching & sorting behavioral frames by trial')

wait_h = waitbar(0,'matching & sorting behavioral frames by trial');
hw=findobj(wait_h,'Type','Patch');
set(hw,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]);

for aa = 1:length(datastruct.pcd.tstamp_mstack)

    
    datastruct.pcd.trial_id(aa) = datastruct.raw.trial_id(idx_vec(aa));
    datastruct.pcd.trial_frame(aa) = datastruct.raw.trial_frame(idx_vec(aa));
    
    wait_h = waitbar((aa) / length(datastruct.pcd.tstamp_mstack));

  
end
disp('finished matching and sorting')
delete(wait_h)

close all

img_obj = nia_movie();
img_obj.loadFlatMovie(cat(3, mean(datastruct.pcd.raw_mstack, 3), datastruct.pcd.raw_mstack));
nia_motionCompensate(img_obj, 1, 1, 1, 1);
corrected_img = img_obj.exportStruct(); 
datastruct.pcd.corr_mstack = nan(size(datastruct.pcd.raw_mstack));
datastruct.pcd.ref_img = corrected_img.slices(1).channels.image;

for ii = 1:length(corrected_img.slices)
        
   datastruct.pcd.corr_mstack(:,:,ii) = corrected_img.slices(ii).channels.image;
        
end

datastruct.pcd.corr_mstack = datastruct.pcd.corr_mstack(:,:,2:end);
        

cd(save_dir)
save('img_dstruct.mat', 'datastruct', '-v7.3')
     
    
cd(homedir)

end
