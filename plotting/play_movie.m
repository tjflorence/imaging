function play_movie(img_struct, fps, do_capture)
close all

img_struct(:,:,1) = img_struct(:,:,2);

num_frames = size(img_struct, 3);
have_tsamps = isfield(img_struct, 'tstamps');

[cmin, cmax] = determine_caxis(img_struct, .6, .001);

if do_capture == 1
mkdir('frames')
cd('frames')
end

if do_capture == 0
    f1 = figure();
else
    f1 = figure('visible', 'off');
    
end

whitebg('black')
set(f1, 'Color', 'k')
for ii = 1:size(img_struct,3)
   
    imagesc(img_struct(:,:,ii));
    set(gca, 'YDir', 'normal')
    axis equal off
    cmap = kjetsmooth(2^16);
    cmap = [0 0 0; cmap];
    colormap(cmap)
    caxis([cmin cmax])
    
   % c1 = colorbar(...
   %     'FontSize', 15, 'FontWeight', 'bold', 'Location', 'South',...
   %     'XTickLabel', {'<100', '150', '200', '250', '300', '>350'});
   %
  % hL = ylabel(c1,'title');     
   % xlabel(c1, 'Pixel Intensity', 'FontSize', 15)
    
    hold on
    
    text(1, 1, num2str(ii))
    
  %  plot([0 10], [5 5], 'Color', 'w', 'LineWidth', 5)
  %  text(-.2,2, '6 um', 'FontSize', 25)
    
    if do_capture ==0
        drawnow
        pause(1/fps)
    end
    
    hold off
    
    if do_capture == 1
        export_fig(['frame_' num2str(ii, '%03d') '.bmp'], '-bmp', '-r140')
    end
    
end
