function plot_dist_vs_activity_matrix(exp_dir, num_divs, roi_num)
% plots distance vs activity matrix and saves
% 

homedir = pwd;
load('roi_data.mat')

% cd to experiment dir. load first training file to get
% first  hyperparams
cd(exp_dir)

exp_files = dir('*train*');
load(exp_files(1).name);
start_i = expr.settings.startXYT(1);
end_i = expr.settings.max_x;

div_length = (end_i-start_i)/num_divs;
div_mat = nan(length(exp_files), num_divs);

try
for ii = 1:(length(exp_files))

    
    load(exp_files(ii).name)
    
    dendvals = expr.c_trial.data.img.roi_trace_raw(roi_num,:);
    
    % find mask for high variance pixels
   % mstack = expr.c_trial.data.img.raw_mstack(:, :, :);
   % mask = make_hivar_mask(mstack, .01);

   % dendvals = nan(1, length(expr.c_trial.data.img.raw_mstack));
   % for aa = 1:length(dendvals)
   %     currframe = mstack(:,:, aa);
   %     pix_i = currframe(mask==1);
    
   %     dendvals(aa) = mean(pix_i);
   % end

    % collect imaging data for each instance of a fly being in a 
    % specific location
    div_vec = nan(1, num_divs);
    start_p = start_i;
    end_p = start_i+div_length;

    for aa = 1:num_divs
        
        % find frame idx when fly is between given position bins
        xpos_inds = find(expr.c_trial.data.xpos>start_p & expr.c_trial.data.xpos<end_p);
   
        % translate each of these into frame 
        img_inds = [];
        for cc = 1:length(xpos_inds)
       
            img_fnum = find(expr.c_trial.data.img.trial_frame==xpos_inds(cc), 1, 'first');
       
                if ~isempty(img_fnum)
                    img_inds = [img_inds img_fnum];
                end
       
        end
   
   
        div_vec(aa) = nanmean(dendvals(img_inds));
   
        start_p = start_p+div_length;
        end_p = end_p+div_length;
   
    end
    
    div_mat(ii, :) = div_vec/(nanmean(div_vec));
end
catch
    disp('incomplete experiment')
end

%div_mat = div_mat/max(max(div_mat));
whitebg('white');
close all

f1 = figure('Color', 'w', 'units', 'normalized', 'Position', [.1 .1 .9 .9])

s1 = subplot(2,1,1)

trench_img = imread('/Users/florencet/Documents/matlab_root/raycast/pink_noise_7_black_wp_1.png');
top_text = trench_img(1:25, :);

image(top_text)
axis equal off
colormap(gray)
hold on

%cpsot = fill([545 675 675 545], [1 1 25.5 25.5], [59 148 255]/255, 'FaceAlpha', .5, 'EdgeAlpha', 0);


freezeColors();

s2 = subplot(2,1,2);
imagesc(div_mat)
hold on
plot( [expr.settings.cool_dist(1)-expr.settings.startXYT(1)...
        expr.settings.cool_dist(1)-expr.settings.startXYT(1)]/div_length, [.5 16.5],...
    'w:', 'LineWidth', 2)

plot( [expr.settings.cool_dist(2)-expr.settings.startXYT(1)...
    expr.settings.cool_dist(2)-expr.settings.startXYT(1)]/div_length, [.5 16.5],...
    'w:', 'LineWidth', 2)


cmap = [[0 0 0]; jet(16)];
colormap(cmap)
caxis([min(min(div_mat)) max(max(div_mat))])
axis equal tight
box off

ticks_here = (expr.settings.max_x-expr.settings.startXYT(1))/4;
tick_vec = [ticks_here 2*ticks_here 3*ticks_here 4*ticks_here]/div_length;
tickmark_vec = ceil(tick_vec*div_length/expr.settings.ticks_per_mm/10);

set(gca, 'XTick', tick_vec, 'XTickLabel', {num2str(tickmark_vec(1)) ...
    num2str(tickmark_vec(2)), num2str(tickmark_vec(3)), num2str(tickmark_vec(4)) });

set(gca, 'YTick', [1 5 10 15 16], 'YTickLabel',{'1' '5' '10', '15', 'probe   '}, 'FontSize', 25)

xlabel('  distance (cm)   ')
ylabel('trial')

c1 = colorbar();
cbar_ticks = linspace(min(min(div_mat)), max(max(div_mat)), 4);
set(c1, 'YTick', cbar_ticks, 'YTickLabel',{'no data   ', ...
    [num2str(cbar_ticks(2)) '   ' ],  [num2str(cbar_ticks(3)) '   ' ],...
    [num2str(cbar_ticks(4)) '   ' ]}, 'FontSize', 25);
ylabel(c1, '  F/bF  ', 'Rotation', -90);


hold on

s2_p = get(s2, 'Position');
s2_pm = [s2_p(1)-.21 s2_p(2)+.25 s2_p(3)*1.61 s2_p(4)*1.1];
set(s2, 'Position', s2_pm)

cd(exp_dir)
mkdir('plots')
cd('plots')
export_fig(['distance_vs_activity_matrix_roi_' num2str(roi_num) '.pdf'], '-pdf', '-zbuffer')

cd(homedir)

end
