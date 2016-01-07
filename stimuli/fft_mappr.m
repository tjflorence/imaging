clear all

blink_fact = 0.5; %% factor in seconds to multiply primes
plist = primes(2); %% list of primes
num_primes = length(plist); %% number of primes / also number of bars

max_cycle_per = blink_fact*plist(end); %% max period of blink cycle
system_hz = 50; %% imaging system frequency
min_num_cycles = 6;

%% id frames to flicker 
command_mat = zeros(length(plist), max_cycle_per*system_hz*min_num_cycles);
for ii = 1:size(command_mat, 1)
    
    c_halfcycle = round(plist(ii)*.5*blink_fact*50);
    mat_p = 0;
    
    while ((mat_p)*c_halfcycle) < size(command_mat, 2)+1
        
        pix_id = (-1)^(mat_p);
        
        if pix_id < 0
            insert_val = 0;
        else
            insert_val = 1;
        end
        
        command_mat(ii,1+(mat_p*c_halfcycle):...
            ((mat_p+1)*c_halfcycle)) = insert_val;
        
        mat_p = mat_p+1;
        
    end
    
    
end

%% set column ids
col_id = nan(32,56);
for ii = 5:5:55
   
    col_id(:,1+(ii-5):ii) = ii/5;
    
end

col_id = circshift(col_id, [1 1]);
num_legal_columns = 11;


%% select columns to populate
randinds = randperm(num_legal_columns);
selected_cols = randinds(1:num_primes);
%stored_cols = [0 0];
%while min(abs(diff(selected_cols)))<2 || ismember(selected_cols, stored_cols, 'rows')
%    randinds = randperm(num_legal_columns);
%    selected_cols = randinds(1:num_primes);
%end
%stored_cols = [stored_cols; selected_cols];



c_cols = nan(32,56);
for ii = 1:length(selected_cols)
   
    c_cols(find(col_id==selected_cols(ii))) = ii;
    
end

blank_frames = 4*ones(32,56, size(command_mat,2));

for ii = 1:size(command_mat,2);
   
    c_frame = 4*ones(32,56);
    
    for jj = 1:num_primes
       
        c_frame(find(c_cols==jj)) = 7*command_mat(jj,ii);
        
    end
    blank_frames(:,:,ii) = c_frame;
    
end


%% collect all this info into the structure, 
%% for vertical stripes only
fmappr.settings.primes = plist;
fmappr.settings.blink_fact = blink_fact;
fmappr.settings.max_cycle_period = max_cycle_per;
fmappr.settings.min_num_cycles = 4;
fmappr.settings.system_hz = 50;
fmappr.settings.command_mat = command_mat;
fmappr.settings.column_id = col_id;

fmappr.vstim(1).frames = blank_frames;
fmappr.vstim(1).current_columns = c_cols;

stim_ind = 1;

%colperms = [nchoosek(1:num_legal_columns, 2); fliplr(nchoosek(1:num_legal_columns,2))];
colperms = [nchoosek(1:num_legal_columns, 1)];... fliplr(nchoosek(1:num_legal_columns,2))];
coldist = abs(diff(colperms, [], 2));
%legalperms = colperms(coldist>1,:);
legalperms = colperms;

rand_ind = randperm(size(legalperms,1));

for aa = 1:length(rand_ind)
   
    
    c_ind = rand_ind(aa);
    selected_cols = legalperms(c_ind,:);

    c_cols = nan(32,56);
    for ii = 1:length(selected_cols)
   
        c_cols(find(col_id==selected_cols(ii))) = ii;
    
    end

    blank_frames = 4*ones(32,56, size(command_mat,2));

    for ii = 1:size(command_mat,2);
   
        c_frame = 4*ones(32,56);
    
        for jj = 1:num_primes
       
            c_frame(find(c_cols==jj)) = 7*command_mat(jj,ii);
        
        end
        blank_frames(:,:,ii) = c_frame;
    
    end

    fmappr.vstim(aa).frames = blank_frames;
    fmappr.vstim(aa).current_columns = c_cols;
    
end

num_legal_rows = 6;
row_id = nan(32, 56);
for ii = 5:5:30

   row_id(ii-3:ii+1,:) = ii/5;

end
stored_rows = [0 0];


%rowperms = [nchoosek(1:num_legal_rows, 2); fliplr(nchoosek(1:num_legal_rows,2))];
rowperms = [nchoosek(1:num_legal_rows, 1)];... fliplr(nchoosek(1:num_legal_rows,2))];
rowdist = abs(diff(rowperms, [], 2));
%legalperms = rowperms(rowdist>1,:);
legalperms = colperms;

rand_ind = randperm(size(legalperms,1));

%% now repeat the process for horizontal stripes
for aa = 1:length(rand_ind)


    c_ind = rand_ind(aa);
    selected_rows = legalperms(c_ind,:);
    
    c_rows = nan(32,56);
    for ii = 1:length(selected_rows)

        c_rows(find(row_id==selected_rows(ii))) = ii;

    end

    blank_frames = 4*ones(32,56, size(command_mat,2));

    for ii = 1:size(command_mat,2);

        c_frame = 4*ones(32,56);
        for jj = 1:num_primes

            c_frame(find(c_rows==jj)) = 7*command_mat(jj,ii);

        end
        blank_frames(:,:,ii) = c_frame;
    end

    fmappr.hstim(aa).frames = blank_frames;
    fmappr.hstim(aa).current_rows = c_rows;
end

cd('/Users/florencet/Documents/matlab_root/imaging/stimuli/')
mkdir('fmappr_frames')
cd('fmappr_frames')

f1 = figure('visible', 'off');
framenum = 1;
for aa = 1:10
   
    for ii = 1:5:size(fmappr.vstim(aa).frames,3)
       
        imagesc(fmappr.vstim(aa).frames(:,:,ii))
        axis equal off
        colormap(gray)
        caxis([0 7])
        
        export_fig(['frame_' num2str(framenum, '%09d') '.bmp'], '-nocrop', '-zbuffer')
        framenum = framenum+1;
        
    end
    
    for ii = 1:5:size(fmappr.hstim(aa).frames,3)
       
        imagesc(fmappr.hstim(aa).frames(:,:,ii))
        axis equal off
        colormap(gray)
        caxis([0 7])
        
        export_fig(['frame_' num2str(framenum, '%09d') '.bmp'], '-nocrop', '-zbuffer')
        framenum = framenum+1;
        
    end
end



v_vec = ones(1, length(blank_frames));
h_vec = 2*ones(1, length(blank_frames));

v_inds = ones(1,length(fmappr.vstim));
h_inds = 2*ones(1, length(fmappr.hstim));
sum_inds = [v_inds, h_inds];
rand_ind = randperm(length(sum_inds));

rndm_vec = sum_inds(rand_ind);

vstim_id = 1;
hstim_id = 1;

type_summary = [];
id_summary = [];
frame_summary = [];

for aa = 1:length(rndm_vec)

    stim_type = rndm_vec(aa);
    type_summary = [type_summary, stim_type*ones(1,length(blank_frames))];
    
    if stim_type == 1
        frame_summary = cat(3, frame_summary, fmappr.vstim(vstim_id).frames);
        id_summary = [id_summary, vstim_id*ones(1, length(blank_frames))];
        vstim_id = vstim_id+1;
        
    elseif stim_type == 2
        
        frame_summary = cat(3, frame_summary, fmappr.hstim(hstim_id).frames);
        id_summary = [id_summary, hstim_id*ones(1, length(blank_frames))];
        hstim_id = hstim_id+1;
        
    end
    
end

fmappr.frame_summary = frame_summary;
fmappr.type_summary = type_summary;
fmappr.id_summary = id_summary;

cd('/Users/florencet/Documents/matlab_root/imaging/stimuli/')
save('fmappr_frames.mat', 'fmappr')
