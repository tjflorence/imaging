function datastruct = exp2struct(varargin)

%% varargs
% arg 1 = imaging directory
% arg 2 = hdf5 directory
% arg 3 = save dir
expdir = varargin{1};
syncdir = varargin{2};


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

frame_out = h5read(h5_f, '/DI/Frame Out');

frame_acq = diff(double(frame_out));
is_frame = find(frame_acq==2);
iframe_int = diff(is_frame);

% have to fudge - missing frame_out signals 
%frame_tstamps = double(is_frame)./acq_hz;
frame_tstamps = linspace(0, length(trial)/acq_hz, num_frames);

%% make frame ministacks
% here i'm assuming 8 frames / stack with 1 extra frame at start

num_zslice = num2str(xml_s.ThorImageExperiment.ZStage.Attributes.steps);
stack_i = 1;
stack_n = num_zslice+1;
data_p = 1;

mstack_mat = nan(frame_x, frame_y, ceil(num_frames/stack_n));
%mstack_mat = nan(frame_x/2, frame_y/2, ceil(num_frames/stack_n));
tstamp_vec = nan(1, ceil(num_frames/stack_n));

wait_h = waitbar(0,'generating ministacks');

while stack_n < num_frames && stack_n < numel(frame_tstamps)
    
    meanframe = sum(frame_mat(:,:,stack_i:stack_n), 3);
  %  clean_frame = medfilt2(meanframe, [8 8]);
    
    ds_frame = meanframe;

    %ds_frame = imresize(meanframe, .5, 'bilinear');
    tstamp_n = frame_tstamps(stack_n);
    
    mstack_mat(:,:,data_p) = ds_frame;
    tstamp_vec(data_p) = tstamp_n;
    
    stack_i = stack_n+1;
    stack_n = stack_i+num_zslice;
    
    data_p = data_p+1;
    
    waitbar( (stack_n) / (num_frames))
    
end

close(wait_h)


%% collect everything into a neat structure
datastruct.pcd.raw_mstack = mstack_mat;
clear mstack_mat

datastruct.pcd.tstamp_mstack = tstamp_vec;
clear tstamp_vec

last_frame = find(isnan(datastruct.pcd.tstamp_mstack)==1, 1, 'first');
datastruct.pcd.tstamp_mstack = datastruct.pcd.tstamp_mstack(1:last_frame);
datastruct.pcd.raw_mstack = datastruct.pcd.raw_mstack(:,:,1:last_frame);

datastruct.raw.pz_pos   = pz_pos(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.clk      = clk(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.trial    = trial(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.frame_out = frame_out(1:(frame_tstamps(end)*acq_hz));
datastruct.raw.timer_stamp  = (1:length(datastruct.raw.pz_pos))/acq_hz;

datastruct.raw.timer_stamp = round2(datastruct.raw.timer_stamp, .0001);
datastruct.pcd.tstamp_mstack = round2(datastruct.pcd.tstamp_mstack, .0001);

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
    
   % diff_vec = abs(datastruct.raw.timer_stamp-stamp);
    %raw_ind = find(diff_vec==min(diff_vec), 1, 'first');
    
    raw_ind = find(datastruct.raw.timer_stamp==stamp);
    
    datastruct.pcd.trial_id(aa) = datastruct.raw.trial_id(raw_ind);
    datastruct.pcd.trial_frame(aa) = datastruct.raw.trial_frame(raw_ind);
    
end



if nargin == 3
   
    cd(varargin{3})
    save('img_dstruct.mat', 'datastruct')
    
end
    
    
cd(homedir)

