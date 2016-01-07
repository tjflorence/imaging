function parse_verification

load('img_dstruct.mat')
homedir = pwd;

whitebg('black');
close all

acq_hz = 10000;

f1 = figure('units', 'normalized', 'position', [.1 .1 .5 .9]);

%% first two subplots are checks to make sure ministacks are being generated
%% properly

% first subplot: number of frames in ministack
s1 = subplot(4,3,1);
hist(datastruct.pcd.verify_vec(:,3)', 100,'b');
box off
xlabel('# frames / ministack', 'fontsize', 15)
ylabel('count', 'fontsize', 15)

s1_p = get(s1, 'Position');
set(s1, 'Position', [s1_p(1)+.08 s1_p(2) s1_p(3) s1_p(4)])

% second subplot: interstack interval
s2 = subplot(4,3,3);
hist(diff(datastruct.pcd.tstamp_mstack), 100)
xlabel('interstack interval (sec)', 'fontsize', 15)

s2_p = get(s2, 'Position');
set(s2, 'Position', [s2_p(1)-.08 s2_p(2) s2_p(3) s2_p(4)])

%% subplots 3-5 are waveforms of piezo position, and frames
total_samples = length(datastruct.raw.pz_pos);
% find random second starts
rand_sample_inds = randperm(total_samples);
choose_seconds = rand_sample_inds(1:3);
sample_separation = diff([choose_seconds total_samples]);
ismore_10000frames = sample_separation> 10000;

% if not in increasing order, separated by > 1 second, and >1 second from
% end, re-roll
while sum(ismore_10000frames) < 3
    
    rand_sample_inds = randperm(total_samples);
    choose_seconds = rand_sample_inds(1:3);
    sample_separation = diff([choose_seconds total_samples]);
    ismore_10000frames = sample_separation> 10000;

end

chosen_seconds = choose_seconds/acq_hz;
frame_acq = diff(datastruct.raw.frame_out);
frame_acq(frame_acq<0) = 0;


% third subplot: random second 1;
frames_of_interest = frame_acq(choose_seconds(1):choose_seconds(1)+acq_hz/2);
pz_of_interest = datastruct.raw.pz_pos(choose_seconds(1):choose_seconds(1)+acq_hz/2);
frame_out_of_interest = datastruct.raw.frame_out(choose_seconds(1):choose_seconds(1)+acq_hz/2);
tstamp_of_interest = datastruct.raw.timer_stamp(choose_seconds(1):choose_seconds(1)+acq_hz/2);

stackinds_of_interest = find(datastruct.pcd.verify_vec(:,1)>choose_seconds(1)...
    & datastruct.pcd.verify_vec(:,2)<choose_seconds(1)+acq_hz/2);

startIds = datastruct.pcd.verify_vec(:,1);
endIds = datastruct.pcd.verify_vec(:,2);
stackind_start_time = datastruct.raw.timer_stamp(startIds(stackinds_of_interest));
stackind_end_time = datastruct.raw.timer_stamp(endIds(stackinds_of_interest));

s3_p = subplot(4,3,4:6);

trace_colors = [255, 67, 250;...
                93, 255, 233]./255;

plot(tstamp_of_interest, frames_of_interest, 'w', 'linewidth', .5)   
hold on

plot(tstamp_of_interest, pz_of_interest, 'color', [.7 .7 .7])

trace_colors = [255, 67, 250;...
                93, 255, 233]./255;

for ii = 1:length(stackinds_of_interest)

    t1 = datastruct.raw.timer_stamp(datastruct.pcd.verify_vec(stackinds_of_interest(ii),1):...
        datastruct.pcd.verify_vec(stackinds_of_interest(ii),2));
    pz_1 = datastruct.raw.pz_pos(datastruct.pcd.verify_vec(stackinds_of_interest(ii),1):...
        datastruct.pcd.verify_vec(stackinds_of_interest(ii),2));
    
    cval = (-1)^ii;
    if cval < 1
        to_plot = trace_colors(1,:);
    else
        to_plot = trace_colors(2,:);
    end
    
    plot(t1, pz_1, 'color', to_plot, 'linewidth', 2)
    
end

xlim([tstamp_of_interest(1) tstamp_of_interest(end)])
ylim([ .99*min(pz_of_interest) 1.01*max(pz_of_interest)]);

box off

% fourth subplot: another random half-second
frames_of_interest = frame_acq(choose_seconds(2):choose_seconds(2)+acq_hz/2);
pz_of_interest = datastruct.raw.pz_pos(choose_seconds(2):choose_seconds(2)+acq_hz/2);
frame_out_of_interest = datastruct.raw.frame_out(choose_seconds(2):choose_seconds(2)+acq_hz/2);
tstamp_of_interest = datastruct.raw.timer_stamp(choose_seconds(2):choose_seconds(2)+acq_hz/2);

stackinds_of_interest = find(datastruct.pcd.verify_vec(:,1)>choose_seconds(2)...
    & datastruct.pcd.verify_vec(:,2)<choose_seconds(2)+acq_hz/2);

startIds = datastruct.pcd.verify_vec(:,1);
endIds = datastruct.pcd.verify_vec(:,2);
stackind_start_time = datastruct.raw.timer_stamp(startIds(stackinds_of_interest));
stackind_end_time = datastruct.raw.timer_stamp(endIds(stackinds_of_interest));

s3_p = subplot(4,3,7:9);

trace_colors = [255, 67, 250;...
                93, 255, 233]./255;

plot(tstamp_of_interest, frames_of_interest, 'w', 'linewidth', .5)   
hold on

plot(tstamp_of_interest, pz_of_interest, 'color', [.7 .7 .7])

trace_colors = [255, 67, 250;...
                93, 255, 233]./255;

for ii = 1:length(stackinds_of_interest)

    t1 = datastruct.raw.timer_stamp(datastruct.pcd.verify_vec(stackinds_of_interest(ii),1):...
        datastruct.pcd.verify_vec(stackinds_of_interest(ii),2));
    pz_1 = datastruct.raw.pz_pos(datastruct.pcd.verify_vec(stackinds_of_interest(ii),1):...
        datastruct.pcd.verify_vec(stackinds_of_interest(ii),2));
    
    cval = (-1)^ii;
    if cval < 1
        to_plot = trace_colors(1,:);
    else
        to_plot = trace_colors(2,:);
    end
    
    plot(t1, pz_1, 'color', to_plot, 'linewidth', 2)
    
end

xlim([tstamp_of_interest(1) tstamp_of_interest(end)])
ylim([ .99*min(pz_of_interest) 1.01*max(pz_of_interest)]);

box off

ylabel('random half-second epochs', 'fontsize', 30)

% fifth subplot: another random half-second
frames_of_interest = frame_acq(choose_seconds(3):choose_seconds(3)+acq_hz/2);
pz_of_interest = datastruct.raw.pz_pos(choose_seconds(3):choose_seconds(3)+acq_hz/2);
frame_out_of_interest = datastruct.raw.frame_out(choose_seconds(3):choose_seconds(3)+acq_hz/2);
tstamp_of_interest = datastruct.raw.timer_stamp(choose_seconds(3):choose_seconds(3)+acq_hz/2);

stackinds_of_interest = find(datastruct.pcd.verify_vec(:,1)>choose_seconds(3)...
    & datastruct.pcd.verify_vec(:,2)<choose_seconds(3)+acq_hz/2);

startIds = datastruct.pcd.verify_vec(:,1);
endIds = datastruct.pcd.verify_vec(:,2);
stackind_start_time = datastruct.raw.timer_stamp(startIds(stackinds_of_interest));
stackind_end_time = datastruct.raw.timer_stamp(endIds(stackinds_of_interest));

s3_p = subplot(4,3,10:12);

trace_colors = [255, 67, 250;...
                93, 255, 233]./255;

plot(tstamp_of_interest, frames_of_interest, 'w', 'linewidth', .5)   
hold on

plot(tstamp_of_interest, pz_of_interest, 'color', [.7 .7 .7])

trace_colors = [255, 67, 250;...
                93, 255, 233]./255;

for ii = 1:length(stackinds_of_interest)

    t1 = datastruct.raw.timer_stamp(datastruct.pcd.verify_vec(stackinds_of_interest(ii),1):...
        datastruct.pcd.verify_vec(stackinds_of_interest(ii),2));
    pz_1 = datastruct.raw.pz_pos(datastruct.pcd.verify_vec(stackinds_of_interest(ii),1):...
        datastruct.pcd.verify_vec(stackinds_of_interest(ii),2));
    
    cval = (-1)^ii;
    if cval < 1
        to_plot = trace_colors(1,:);
    else
        to_plot = trace_colors(2,:);
    end
    
    plot(t1, pz_1, 'color', to_plot, 'linewidth', 2)
    
end

xlim([tstamp_of_interest(1) tstamp_of_interest(end)])
ylim([ .99*min(pz_of_interest) 1.01*max(pz_of_interest)]);

box off
xlabel('time (sec)', 'FontSize', (30))

mkdir('plots')
cd('plots')

export_fig('ministack_parsing.pdf', '-pdf', '-zbuffer', '-r140')
cd(homedir)














