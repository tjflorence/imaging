function roi_select(trial_num)

home_dir = pwd;

exp_trials = dir('env*');
load(exp_trials(1).name)

cMap = [255, 70, 69; ...
        186, 48, 232; ...
        65, 76, 255; ...
        48, 201, 232; ...
        24, 255, 103; ...
        232, 232, 45; ...
        232, 131, 17]./255;
    
cMap = repmat(cMap, [100, 1]);

neuron_fig = figure();
whitebg('black')
imagesc(std(expr.c_trial.data.img.corr_mstack,[], 3));
caxis([3 20])
colormap([0,0,0; 0,0,0; jet(8)])
axis equal off

roi_response = 1;
ii = 1;

while roi_response ~= 0
  
  mask_response = 0;
  while mask_response ~= 1
    
    disp('select next ROI');
    [BW, xi, yi] = roipoly;
  
    roi_struct(ii).xy = [xi, yi];
  

    
    temp_mask_fig = figure();
    imagesc(BW)
    colormap(gray)
    axis equal off
  
    mask_response = input('accept mask? 1 = yes, 2 = no');
    close(temp_mask_fig)
  
  end
  
  roi_struct(ii).mask = BW;
  roi_struct(ii).cmap = cMap(ii,:);
  hold on
  
  scat_h = fill(xi, yi, cMap(ii,:));
  set(scat_h, 'MarkerFaceColor', cMap(ii,:), 'MarkerEdgeColor', 'w');
  
  roi_response = input('additional ROI? 1 = yes, 0 = no (complete)');
  
  ii = ii+1;
  
    
    
end

save('roi_data.mat', 'roi_struct')

mkdir('plots')
cd('plots')
close all

neuron_fig = figure();
whitebg('black')
imagesc(mean(expr.c_trial.data.img.corr_mstack, 3));
colormap([0,0,0; 0,0,0; gray(8)])
axis equal off

hold on
for jj = 1:length(roi_struct)
    
  scat_h = fill(roi_struct(jj).xy(:,1), roi_struct(jj).xy(:,2), cMap(jj,:));
  set(scat_h, 'MarkerFaceColor', cMap(jj,:), 'MarkerEdgeColor', 'w');
  
end

export_fig('selected_roi.pdf', '-pdf', '-zbuffer')

cd(home_dir);
end

