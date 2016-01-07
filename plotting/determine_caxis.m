function  [cmin, cmax] = determine_caxis(img_data, low_pct, high_pct)


frame1 = img_data(:,:,1);
framepix_vec = frame1(:);

[n, xout] = hist(framepix_vec, 1000);

n_flip = fliplr(n);
xout_flip = fliplr(xout);


norm_n = cumsum(n/sum(n));
flip_norm_n = cumsum(n_flip/sum(n));


low_idx = find(norm_n>low_pct, 1, 'first');
hi_idx = find(flip_norm_n> high_pct, 1, 'first');

cmin = xout(low_idx);
cmax = xout_flip(hi_idx);

end
