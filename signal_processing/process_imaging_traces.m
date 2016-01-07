function process_imaging_traces(expfile, crop_i, crop_n)
% PROCESS_IMAGING_TRACES grabs single trace from imaging data. 
% EXPFILE is expr file after initial processing
% CROP_I is start index of crop on imaging data
% CROP_N is end index of crop on imaging data
% saves updated file to disk

    load(expfile)

    % pull out mstack, crop to desired size
    mstack = expr.c_trial.data.img.raw_mstack;
    mstack = mstack(crop_i:crop_n, crop_i:crop_n, :);

    % get high variance mask
    mask = make_hivar_mask(mstack, .05);
    num_frames = size(mstack,3);
    raw_frameval_vec = nan(1, num_frames);

    % for each frame, get mean value of high-variance pixels
    for ii = 1:num_frames
    
        c_frame = mstack(:,:,ii);
        frame_vals = mean(c_frame(mask==1));
    
        raw_frameval_vec(ii) = frame_vals;

    end

    % find when the lights turn on at each trial; dark frames will be our
    % baseline
    frame_at_dark = expr.c_trial.data.frame_p(expr.c_trial.data.state_1_2_trans);
    df_trace = (raw_frameval_vec-mean(raw_frameval_vec(1:frame_at_dark))) /...
                mean(raw_frameval_vec(1:frame_at_dark)) ;

    expr.c_trial.data.img.raw_trace = raw_frameval_vec;
    expr.c_trial.data.img.df_trace = df_trace;
    
    % save it out
    save(expfile, 'expr')

end