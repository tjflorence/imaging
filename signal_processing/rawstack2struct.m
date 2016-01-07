function rawstack2struct(expdir, syncdir, savedir)
%% varargs
% arg 1 = imaging directory
% arg 2 = hdf5 directory
% arg 3 = save dir
% description: uses piezo position to make experiment ministacks, saves as
% structure

acq_hz = 10000;
homedir = pwd;

%% load thorimage imaging data - binary file and xml parameters
%expdir = varargin{1};
cd(expdir)

img_f = 'Image_0001_0001.raw';
xml_f = 'Experiment.xml';

xml_s = xml2struct(xml_f);
frame_x = str2num(xml_s.ThorImageExperiment.LSM.Attributes.pixelX);
frame_y = str2num(xml_s.ThorImageExperiment.LSM.Attributes.pixelY);
disp(['frame dimensions: ' num2str(frame_y) ' x ' num2str(frame_x)])

frame_dim = [frame_x frame_y];

disp('loading RAW file')
f_h = fopen(img_f, 'r');
raw_out = fread(f_h, 'uint16');
fclose(f_h);
disp('RAW file loaded')

num_frames = double(size(raw_out,1)/(frame_dim(1)*frame_dim(2)));
frame_mat = reshape(raw_out, [frame_dim(1), frame_dim(2), num_frames]);

%% load thorsync
%syncdir = varargin{2};
cd(syncdir)
h5_f = 'Episode001.h5';

pz_pos  = h5read(h5_f, '/AI/Piezo Monitor');
clk     = h5read(h5_f, '/AI/clk');
trial   = h5read(h5_f, '/AI/trial');
try
    measured_temp = h5read(h5_f, '/AI/set_temp'); % these are flipped, need to fix
    set_temp = h5read(h5_f, '/AI/measured_temp');
catch
    disp('oops, no temperature data')
end
frame_out = h5read(h5_f, '/DI/Frame Out');

% differentiate frame_out to find frame_acq frames
frame_acq = diff(double(frame_out));
is_frame = find(frame_acq>1.9);

% take frame_acq and make frame_counter - this will allow us to go from 
% index to frame count later
disp('generating frame ids')
fcounter = 0;
f_countid = zeros(1,length(frame_out));
for ii = 1:length(f_countid)-1
    if frame_acq(ii)>1.9
        fcounter = fcounter+1;
    end
    f_countid(ii) = fcounter;
end
f_countid(end) = f_countid(end-1);

if numel(is_frame) ~= num_frames
    disp('# sync frames does NOT equal # stacks')
else
    disp('# sync frames == # stacks')
end

% have to fudge - missing frame_out signals 
%frame_tstamps = double(is_frame)./acq_hz;
frame_tstamps = linspace(0, length(trial)/acq_hz, num_frames);
sync_tstamps = 0:(1/acq_hz):(length(frame_out)/acq_hz);

%% make frame ministacks, using pz position
% step 1: smooth pz position by convolving with gaussian
% we are acquiring signals at 10Khz - so smooth with window of 300 samples?
% step 2: invert and find peaks
w = 300;
sigma = 10;
x = linspace(-w / 2, w / 2, w);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); 
pz_sm = conv(pz_pos, gaussFilter, 'same');
% we're filtered - invert waveform and find peaks
pz_sm = -pz_sm;
[pks, loc] = findpeaks(pz_sm);

% still get spurious peaks. filter these out
[count, val] = hist(pks);
count(count==0) = nan;
min_pos = val(find(count==min(count)));
max_pos = val(find(count==max(count)));
mid_pos = mean([min_pos, max_pos]);

filt_pks = pks(pks>mid_pos);
filt_loc = loc(pks>mid_pos);

% filt_loc is the filtered indices of piezo troughs
% all frames within filt_loc indices are members of same z_stack

num_zslice = num2str(xml_s.ThorImageExperiment.ZStage.Attributes.steps);

% pre-allocate matrices
mstack_mat = nan(frame_x, frame_y, ceil(num_frames/num_zslice));
tstamp_vec = nan(1, ceil(num_frames/num_zslice));
mtemp_vec = nan(1, ceil(num_frames/num_zslice));
stemp_vec = nan(1, ceil(num_frames/num_zslice));
verify_vec = nan(ceil(num_frames)/num_zslice, 3); % verify_vec exists to keep trace of which frames ended up in which stack.
% verify_vec = [ministack_start ministack_end num_slices]

wait_h = waitbar(0,'generating ministacks');

stack_i = 1;
stack_n = filt_loc(1);
data_p = 0;

grab_frames = 0;
stoploop = 0;
last_time = 0;

for ii = 1:length(filt_loc)
    
    stack_n = filt_loc(ii);
    grab_frames = f_countid(stack_i:stack_n);
    
    if grab_frames(1) == 0
        grab_frames(1) = 1;
    end
    
    tstamp_n = sync_tstamps(stack_n);
    try
        mtemp_n = measured_temp(stack_n);
        stemp_n = set_temp(stack_n);
    catch
        mtemp_n = nan;
        stemp_n = nan;
    end

    if (tstamp_n - last_time) > .005
    
    
        meanframe = sum(frame_mat(:,:,grab_frames(1):grab_frames(end)), 3);    
        ds_frame = meanframe;

        data_p = data_p+1;
    
        mstack_mat(:,:,data_p) = ds_frame;
        tstamp_vec(:,data_p) = tstamp_n;
        mtemp_vec(:,data_p) = mtemp_n;
        stemp_vec(:,data_p) = stemp_n;
        verify_vec(data_p,:) = [stack_i stack_n (grab_frames(end)-grab_frames(1))]; 
    
        stack_i = stack_n;
        last_time = tstamp_n;

    end
    
    wait_h = waitbar( (ii) / length(filt_loc));
    
end

close(wait_h)
close all

%% collect everything into a neat structure
datastruct.pcd.raw_mstack = mstack_mat;
datastruct.pcd.tstamp_mstack = tstamp_vec;
datastruct.pcd.measured_temps = mtemp_vec;
datastruct.pcd.set_temps = stemp_vec;
datastruct.pcd.verify_vec = verify_vec;
datastruct.pcd.verify_vec_id = {'miniframe_start', 'miniframe_end', 'num_frames'};
clear mstack_mat mtemp_vec stemp_vec verify_vec tstamp_vec


datastruct.raw.pz_pos   = pz_pos(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.clk      = clk(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.trial    = trial(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.frame_out = frame_out(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.timer_stamp  = sync_tstamps;
datastruct.raw.measured_temps = measured_temp;
datastruct.raw.set_temps = set_temp;

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

for ii = 1:length(t_ends)
   
    t_start_i = t_starts(ii);
    t_end_i = t_ends(ii);
    
    datastruct.raw.trial_id(t_start_i:t_end_i) = ii;
    
end


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

end

datastruct.pcd.trial_id = nan(1, length(datastruct.pcd.tstamp_mstack));
datastruct.pcd.frame_id = nan(1, length(datastruct.pcd.tstamp_mstack));

for aa = 1:length(datastruct.pcd.tstamp_mstack)-1
    stamp = datastruct.pcd.tstamp_mstack(aa);
    
    %diff_vec = abs(datastruct.raw.timer_stamp-stamp);
    %raw_ind = find(diff_vec==min(diff_vec), 1, 'first');
    
    raw_ind = find(datastruct.raw.timer_stamp==stamp);
    
    datastruct.pcd.trial_id(aa) = datastruct.raw.trial_id(raw_ind);
    datastruct.pcd.trial_frame(aa) = datastruct.raw.trial_frame(raw_ind);
    
end



cd(savedir)
save('img_dstruct.mat', 'datastruct', '-v7.3')
     
    
cd(homedir)

