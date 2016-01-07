function img2behav(expdir)
%% new comment

num_reframes = 10;

homedir = pwd;

cd(expdir)

cd(expdir)
load('img_dstruct.mat')

bfiles = dir('env_*');

for aa = 1:length(bfiles)

    load(bfiles(aa).name)
    disp(['running trial ' num2str(aa) ' out of ' num2str(length(bfiles))])
    
    exp_frames = find(datastruct.pcd.trial_id==aa);
    expr.c_trial.data.img.raw_mstack = datastruct.pcd.raw_mstack(:,:,exp_frames);
    expr.c_trial.data.img.corr_mstack = datastruct.pcd.corr_mstack(:,:,exp_frames);

    expr.c_trial.data.img.trial_frame = datastruct.pcd.trial_frame(exp_frames);
    expr.c_trial.data.img.mtemps = datastruct.pcd.measured_temps(exp_frames);
    expr.c_trial.data.img.stemps = datastruct.pcd.set_temps(exp_frames);
    

    expr.c_trial.data.img.mean_img = mean(datastruct.pcd.raw_mstack(:,:,1), 3);
    expr.c_trial.data.img.ref_count = Inf;
      
    
    
    expr.c_trial.data.frame_p = nan(1,expr.c_trial.data.count);
    enum = 1;
    for bb = 1:expr.c_trial.data.count
       
        if bb < expr.c_trial.data.img.trial_frame(enum)
            expr.c_trial.data.frame_p(bb) = expr.c_trial.data.img.trial_frame(enum);
            
        else
            
            while bb > expr.c_trial.data.img.trial_frame(enum) &&...
                    enum < numel(expr.c_trial.data.img.trial_frame)
                enum = enum+1;
            end
            
            expr.c_trial.data.frame_p(bb) = expr.c_trial.data.img.trial_frame(enum);
            
        end
        
          
    end
    
    clear img_obj corrected_img
    
    save(bfiles(aa).name, 'expr', '-v7.3')

end

cd(homedir)

end
