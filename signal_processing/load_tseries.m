function tiff_struct = load_tseries(dir_name)

cd(dir_name)

%% get tif stacks
tif_files = dir('ChanA_0*');
num_frames = length(tif_files);

%% handle to timing file
timing_file = 'timing.txt';
tfile_h = fopen(timing_file);

%% initialize frame matrix, timing matrix
test_frame = imread(tif_files(1).name);
size_mn = size(test_frame);

tif_mat = nan(size_mn(1), size_mn(2), num_frames);
timing_vec = [];

%% read in tifs to matrix
for ii = 1:num_frames
    
    tif_mat(:,:,ii) = imread(tif_files(ii).name);

end

%% read timing file lines to timing vec
out = 0;
% first line is zero
out = fgetl(tfile_h);
while out ~= -1
   
    out = fgetl(tfile_h);
    timing_vec = [timing_vec str2double(out)];
    
end

fclose(tfile_h);

%% check that timing frames, tifs have same number of samples
num_tstamps = length(timing_vec);

if num_frames == num_tstamps
    disp('number of frames and timestamps are equal')
else
    disp('number of frames and timestamps do not match')
end


tiff_struct.img = tif_mat;
tiff_struct.tstamps = timing_vec;

end
