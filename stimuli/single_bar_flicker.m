clear all

system_hz = 50;
half_hz = system_hz/2;
num_flickers = 5;

single_cycle = [ones(1, half_hz), zeros(1, half_hz)];
command_vec = repmat(single_cycle, [1 num_flickers]);

sample_frame = nan(32,56);

col_ids = nan(32,56);
colnum_id = 0;
for ii = 1:2:56 % sets inter-bar distance
    
    colnum_id = colnum_id+1;    
    col_ids(:,ii:ii+1) = colnum_id; % sets bar width
    
end

row_ids = nan(32,56);
rownum_id = 0;
for ii = 1:2:32
    
    rownum_id = rownum_id+1;
    row_ids(ii:ii+1,:) = rownum_id;

end

num_cols = colnum_id;
num_rows = rownum_id;

for ii = 1:num_cols
    
    blank_frames = 3*ones(32,56,length(command_vec));
   
    for jj = 1:length(command_vec)
        
        blank_frame = 3*ones(32,56);
        blank_frame(find(col_ids==ii)) = 7*command_vec(jj);
     
        blank_frames(:,:,jj) = blank_frame;
        
    end
    
    flick_stim.vstim(ii).frames = blank_frames;
    flick_stim.vstim(ii).col_id = ii;
    flick_stimm.vstim(ii).type = 1;
    
end

for ii = 1:num_rows
    
    blank_frames = 3*ones(32,56,length(command_vec));
   
    for jj = 1:length(command_vec)
        
        blank_frame = 3*ones(32,56);
        blank_frame(find(row_ids==ii)) = 7*command_vec(jj);
     
        blank_frames(:,:,jj) = blank_frame;
        
    end
    
    flick_stim.hstim(ii).frames = blank_frames;
    flick_stim.hstim(ii).col_id = ii;
    flick_stim.hstim(ii).type = 2;
    
end

type_vec = [ones(1, num_cols) 2*ones(1,num_rows)];
flick_stim.type_vec = type_vec(randperm(length(type_vec)));
frame_summary = [];
type_summary = [];
spot_summary = [];

v_p = 1;
h_p = 1;
rand_cols = randperm(num_cols);
rand_rows = randperm(num_rows);
for ii = 1:length(flick_stim.type_vec)
    
    if flick_stim.type_vec(ii) == 1
       
        col_id = rand_cols(v_p);
        frame_summary = cat(3,frame_summary, flick_stim.vstim(col_id).frames);
        type_summary = [type_summary, ones(1, length(flick_stim.vstim(col_id).frames))];
        spot_summary = [spot_summary, col_id*ones(1, length(flick_stim.vstim(col_id).frames))];
        v_p = v_p+1;
        
    else
        
        row_id = rand_rows(h_p);
        frame_summary = cat(3,frame_summary, flick_stim.hstim(row_id).frames);
        type_summary = [type_summary, 2*ones(1, length(flick_stim.hstim(row_id).frames))];
        spot_summary = [spot_summary, row_id*ones(1, length(flick_stim.hstim(row_id).frames))];
        h_p = h_p+1;
        
    end

end

flick_stim.frame_summary = cat(3,frame_summary, frame_summary, frame_summary, frame_summary);
flick_stim.type_summary = [type_summary type_summary type_summary type_summary];
flick_stim.spot_summary = [spot_summary spot_summary spot_summary spot_summary];


cd('/Users/florencet/Documents/matlab_root/imaging/stimuli/')
save('small_onebar_frames.mat', 'flick_stim')

mkdir('small_onebar_frames')
cd('small_onebar_frames')
f1 = figure('visible', 'off');
framenum = 1;
for aa = 1:10:length(frame_summary)
    
    imagesc(frame_summary(:,:,aa))
    axis equal off
    colormap(gray)
    caxis([0 7])
        
    export_fig(['frame_' num2str(framenum, '%09d') '.bmp'], '-nocrop', '-zbuffer')
    framenum = framenum+1;
   
end




