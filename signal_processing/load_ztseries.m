function tiff_struct = load_ztseries(idirname, hdf5dirname,...
                                     ppstack, index_vec)

curr_dir = pwd;


%% check that indices are appropriate for stack
if mod(index_vec(2) - index_vec(1), ppstack) ~= (ppstack-1)
    err('frame indices are not appropriate for stack')
end

%% load hdf5 file frame timing acquisition data
disp('loading timestamps')
cd(hdf5dirname)
h5_file = dir('*.h5');
frame_in = h5read(h5_file(1).name, '/DI/Frame Out');
disp('timestamps loaded')

lead_edge_ind   = find(diff(frame_in)==2);
lead_edge_time  = lead_edge_ind*(1/30000);

%% get tif stack fnames
cd(idirname)
tif_files = dir('ChanA_0*');
num_frames = ceil( (index_vec(2)-index_vec(1))/ppstack);

%% initialize frame matrix, timing matrix
test_frame = imread(tif_files(1).name);
size_mn = size(test_frame);

tif_mat     = nan(size_mn(1), size_mn(2), num_frames);
stack_mat   = nan(size_mn(1), size_mn(2), ppstack);
timing_vec = [];

%% read in tifs to matrix, tstamps to vec
disp('loading frames')
ministack_ind = 1;
for ii = index_vec(1):ppstack:index_vec(2)
    
    for jj = 0:(ppstack-1)
       
        stack_mat(:,:,jj+1) = imread(tif_files(ii+jj).name);
        
        if jj == ppstack-1
            
            timing_vec = [timing_vec lead_edge_time(ii+jj)];
            
        end
        
    end
    
    tif_mat(:,:,ministack_ind) = mean(stack_mat,3);
    stack_mat   = nan(size_mn(1), size_mn(2), ppstack);
    ministack_ind = ministack_ind+1;
end
disp('frames loaded')

clear tif_files h5_file
cd(curr_dir)

tiff_struct.img = tif_mat;
tiff_struct.tstamps = timing_vec;
tiff_struct.idirname = idirname;
tiff_struct.hdf5dirname = hdf5dirname;

end
