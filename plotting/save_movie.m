function save_movie(img_struct)
close all


movie_name = '.2 um bead wobble'
num_frames = size(img_struct.img, 3);
have_tstamps = isfield(img_struct, 'tstamps');


mkdir('movie_frames')
cd('movie_frames')

f1 = figure();
whitebg('black')
set(f1, 'Color', 'k')

for ii = 1:num_frames
   
    imagesc(img_struct.img(:,:,ii));
    
    axis equal off
    colormap gray
    caxis([0 2500])
    
    
    text(-50, -30, movie_name,...
        'Color', 'w', 'FontSize', 15,...
        'FontWeight', 'bold')
    
    
    frame_num = sprintf('%03d',ii) ;
    
    
    text(-50, 790, ['frame: ' frame_num '                                          '...
       'time: ' num2str(img_struct.tstamps(ii)) ' sec'],...
        'Color', 'w', 'FontSize', 15,...
        'FontWeight', 'bold')
    
    
    export_fig(['movie_frame_' frame_num '.bmp'], '-bmp', '-nocrop')
    %drawnow
    %pause(1/fps)
    
    
end
