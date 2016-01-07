function add_trace_data

exp_files = dir('env*');
load('roi_data.mat')
load(exp_files(2).name)

%expr.c_trial.data.img.roi_trace_raw = nan(length(roi_struct),...
%                                    size(expr.c_trial.data.img.corr_mstack, 3));
                                

%% add roi for whole experiment
load('img_dstruct.mat')
try
    datastruct.pcd = rmfield(datastruct.pcd, 'roi_trace_dF');
    datastruct.pcd = rmfield(datastruct.pcd, 'roi_trace_raw');
    datastruct.pcd = rmfield(datastruct.pcd, 'roi_baseline');
catch
    disp('no old trace data')
end
    
dastastruct.pcd.roi_trace_raw = nan(length(roi_struct), size(datastruct.pcd.corr_mstack, 3));
dastastruct.pcd.roi_trace_dF = nan(length(roi_struct), size(datastruct.pcd.corr_mstack, 3));
dastastruct.pcd.roi_baseline = nan(length(roi_struct), 3);
    
% grab raw pixels from roi's
for aa = 1:length(roi_struct)
       
    for ii = 1:size(datastruct.pcd.corr_mstack, 3)
           c_frame = datastruct.pcd.corr_mstack(:,:,ii);
           roi_pix = c_frame(roi_struct(aa).mask==1);
           
           datastruct.pcd.roi_trace_raw(aa,ii) = mean(roi_pix);
    end
        
        
end
    
for aa = 1:length(roi_struct)
    datastruct.pcd.roi_baseline(aa) = mean(datastruct.pcd.roi_trace_raw(aa,:));
    datastruct.pcd.roi_trace_dF(aa,:) = (datastruct.pcd.roi_trace_raw(aa,:)-datastruct.pcd.roi_baseline(aa)) /datastruct.pcd.roi_baseline(aa);
end

save('img_dstruct.mat', 'datastruct', '-v7.3')

    

                                
%% add roi for each trial                             
for aa = 1:length(exp_files)
    
    load(exp_files(aa).name)
    try
    expr.c_trial.data.img = rmfield(expr.c_trial.data.img, 'roi_trace_dF');
    expr.c_trial.data.img = rmfield(expr.c_trial.data.img, 'roi_trace_raw');
    expr.c_trial.data.img = rmfield(expr.c_trial.data.img, 'roi_baseline');
    catch
        disp('no old trace data')
    end

    for ii = 1:length(roi_struct)
        for jj = 1:size(expr.c_trial.data.img.corr_mstack,3)
        
            c_frame = expr.c_trial.data.img.corr_mstack(:,:,jj);
            roi_pix = c_frame(roi_struct(ii).mask==1);
        
            expr.c_trial.data.img.roi_trace_raw(ii,jj) = mean(roi_pix);  
        
        end
    end
    
    expr.c_trial.data.img.roi_trace_raw(:,1) = expr.c_trial.data.img.roi_trace_raw(:,2);
    
    expr.c_trial.data.img.roi_trace_dF = nan(size(expr.c_trial.data.img.roi_trace_raw));
    
    % now baseline subtract
    expr.c_trial.data.img.roi_baseline = nan(length(roi_struct), 1);
    for ii = 1:length(roi_struct)
        
        c_trace = expr.c_trial.data.img.roi_trace_raw(ii, :);
        
        last_darkframe = find(expr.c_trial.data.img.trial_frame < ...
                                expr.c_trial.data.state_1_2_trans, 1, 'first');
                            
        if ~isnan(last_darkframe)
            disp('dark baseline!')
            baseline_val = mean(c_trace(1:last_darkframe));
        
        else
            disp('no dark baseline')
            baseline_val = mean(c_trace);
            
            
        end
        
        expr.c_trial.data.img.roi_trace_dF(ii, :) = (c_trace-baseline_val)/baseline_val;
        
        expr.c_trial.data.img.roi_baseline(ii) = baseline_val;
        
        
        
    end
    
    save(exp_files(aa).name, 'expr', '-v7.3')
    
    
end

