function mask = make_hivar_mask(img_stack, pctile)
% make_hivar_mask creates mask for high-variance pixels in an image stack
% IMG_STACK is y * x * num_frames matrix
% PCTILE is desired threshold in decimal (ie top 5% is .05)
% returns MASK, a y * x zeros matrix with ones in place of high-variance pixels 

    % make histogram of pixel values
    std_img = std(img_stack, [] , 3);
    [n, xout] = hist(std_img(:), 100);

    % begin integrating
    sum_q = 0;
    p = 2;
    while sum_q < (1-pctile)
    
        sum_q = sum(n(1:p))/sum(n);
        p = p+1;
    
    end
    
    threshval = xout(p);
    
    mask = zeros(size(std_img));
    mask(std_img>threshval)=1;

end