function nia_playGenericMovie(mov, callback_list)
%NIA_PLAYGENERICMOVIE Play the passed generic movie.
%   nia_playGenericMovie(mov, callback_list) displays the passed movie.
%   The argument mov can be any object that the functions of the
%   callback_list are prepared to handle.
%
%   The callback_list argument provides a set of callbacks used
%   to provide information on a generic movie object.  It must
%   have the following fields:
%
%       getImage - Function handle that accepts the movie object,
%           a position vector, a colormap specification, and returns a 2D
%           matrix image and a time value.
%
%       getChannels - Function handle that accepts the movie object,
%           and returns a [1x(number of channels)] array of channel
%           identifiers. The number of channels must be greater than
%           one.
%
%       getPosRanges - Function handle that accepts the movie object,
%           and returns a [2x(number of position dimensions)] array. The
%           first row contains the start of the allowed range for each
%           position, and the second row contains the end of the allowed
%           range, both the start and end of the range are given
%           inclusively. The number of positions must be greater than
%           one.
%
%       getHistInfo - Function handle that accepts the user data and
%           movie object and returns a structure array of histogram
%           information. The required histogram information is described in
%           detail in the documentation for the function
%           nia_selectColormap().
%
%       processROIs - Function handles that accepts the user data and
%           movie object and a structure array of ROIs.  The structure
%           array of ROI contains the fields 'handle' which is an imroi
%           object, the field 'channel' which is the channel selected for
%           the ROI, the field 'home_pos' which is the home position of
%           the ROI, and the field 'scan_pos' which is the position index
%           to be used to produce the series. It must return a cell array
%           with (num ROI) elements of [2x(num positions in scan_pos)]
%           matrices containing the position trace of the average pixel
%           intensity within the ROI. Note that the returned vectors must
%           not be empty. If no valid slices were found, then this function
%           should return a [0; NaN] value.
%

% TODO:
% - should be able to change color and channel of ROIs and
%   get live update in trace window
% - the context menu for vertices of rois is still screwy

% INTERNAL NOTES:
% - The position sliders are used as the authoritative source of
%   the current position.  The retrieved value should be rounded
%   prior to being used in any way.

% Bugs I can't fix because matlab has issues:
% -  the imrect tool is disabled, can't be fixed because the Select
%    property of the imrect object doesn't work consistently
% -  the backgrounds on button icons aren't quite right, can't be
%    fixed because matlab doesn't allow for alpha channel in images


cfg_default_rate = 10;
cfg_timer_interval = 0.03;
cfg_cmap_fig_xpos = 50;
cfg_cmap_fig_ypos = -50;
cfg_cmap_fig_width = 400;
cfg_cmap_fig_height = 400;
cfg_label_font_size = 12;

% Check arguments
if ~isstruct(callback_list) || length(callback_list) ~= 1
    error 'The argument ''callback_list'' has an invalid type';
end

cb_list_names = {'getImage', 'getChannels', 'getPosRanges', ...
    'getHistInfo', 'processROIs'};

[cb_list_ok, cb_list_msg] = nia_hasValidFieldNames(...
    callback_list, cb_list_names, cb_list_names);
if ~cb_list_ok
    error(cb_list_msg, 'callback_list');
end

% Get position ranges so we can configure sliders
pos_ranges = callback_list.getPosRanges(mov);

% Create the figure
fig = figure;
set(fig, 'MenuBar', 'none');
fig_pos = get(fig, 'Position');
fig_bkg = get(fig, 'Color');

% Create the toolbar
% NOTE: handle visibility is set to off below because MATLAB tutorial
% says to do that.  I have no idea why, matlab has a pretty terrible
% event propagation model.
tbar = uitoolbar(fig);

poly_image = imread('nia_poly_icon.png');
poly_ctrl = uipushtool(tbar, ...
    'CData', poly_image, ...
    'TooltipString', 'Draw ROI polygon', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

ellipse_image = imread('nia_ellipse_icon.png');
ellipse_ctrl = uipushtool(tbar, ...
    'CData', ellipse_image, ...
    'TooltipString', 'Draw ROI ellipse', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

plot_image = imread('nia_plot_icon.png');
plot_ctrl = uipushtool(tbar, ...
    'CData', plot_image, ...
    'TooltipString', 'Plot ROI values', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Enable', 'off');

cmap_image = imread('nia_cmap_icon.png');
cmap_ctrl = uipushtool(tbar, ...
    'CData', cmap_image, ...
    'TooltipString', 'Configure colormaps', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Enable', 'on');

% Create the play button
play_btn_play_image = imread('nia_play_icon.png');
play_btn_pause_image = imread('nia_pause_icon.png');

% Create position sliders
pos_text_labels = zeros(1, size(pos_ranges, 2));
pos_play_btns = zeros(1, size(pos_ranges, 2));
pos_slider_ctrls = zeros(1, size(pos_ranges, 2));
pos_num_edits = zeros(1, size(pos_ranges, 2));

for pos_idx=1:size(pos_ranges,2)
    pos_step = pos_ranges(2, pos_idx) - ...
        pos_ranges(1, pos_idx);
    
    % MATLAB requires that min value of slider is less than
    % max value.  That is probably a bug, but we can work around
    % it by bumping the maximum value to avoid exact equality.
    if pos_step == 0
        pos_step = 1;
        pos_bump = 1e-3;
    else
        pos_bump = 0;
    end
    
    pos_text_labels(pos_idx) = uicontrol(...
        'Style', 'text', ...
        'Units', 'pixels', ...
        'String', sprintf('Position %d', pos_idx), ...
        'FontSize', cfg_label_font_size, ...
        'HorizontalAlignment', 'left', ...
        'Background', fig_bkg);

    pos_play_btns(pos_idx) = uicontrol(...
        'Style', 'pushbutton', ...
        'Units', 'pixels', ...
        'BackgroundColor', fig_bkg, ...
        'CData', play_btn_play_image);
    
    pos_slider_ctrls(pos_idx) = uicontrol(...
        'Style', 'slider', ...
        'Units', 'pixels', ...
        'Value', pos_ranges(1, pos_idx), ...
        'Min', pos_ranges(1, pos_idx), ...
        'Max', pos_ranges(2, pos_idx) + pos_bump, ...
        'SliderStep', [1/pos_step, 0.05]);
    
    pos_num_edits(pos_idx) = uicontrol(...
        'Style', 'edit', ...
        'Units', 'pixels', ...
        'String', num2str(pos_ranges(1, pos_idx)), ...
        'BackgroundColor', fig_bkg, ...
        'Min', 0, ...
        'Max', 1);
end

% Create the framerate edit
rate_label = uicontrol(...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'String', 'Framerate', ...
    'FontSize', cfg_label_font_size, ...
    'HorizontalAlignment', 'left', ...
    'Background', fig_bkg);

rate_ctrl = uicontrol(...
    'Style', 'edit', ...
    'Units', 'pixels', ...
    'String', num2str(cfg_default_rate), ...
    'Background', fig_bkg, ...
    'Min', 0, ...
    'Max', 1);

% Create the timelabel
time_label = uicontrol(...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'String', 'Time: ', ...
    'FontSize', cfg_label_font_size, ...
    'HorizontalAlignment', 'right', ...
    'Background', fig_bkg);

% Create the main image axes
mov_ax = axes;
set(mov_ax, ...
    'Unit', 'pixels', ...
    'XTick', [], ...
    'YTick', [], ...
    'Box', 'on');

% Create the colormap selector
udata.colormap_spec = [];
set(fig, 'UserData', udata);

cmap_fig = figure('Visible', 'off');
set(cmap_fig, 'Position', [fig_pos(1) + cfg_cmap_fig_xpos, ...
    fig_pos(2) + cfg_cmap_fig_ypos, cfg_cmap_fig_width, ...
    cfg_cmap_fig_height]);

ch_ids = callback_list.getChannels(mov);
hist_info = callback_list.getHistInfo(mov);
udata.colormap_spec = nia_selectColormap(...
    cmap_fig, ch_ids, hist_info, ...
    @(spec) colormapChanged(fig, spec));
set(fig, 'UserData', udata);

% Create the playback timer
timer_obj = timer();
set(timer_obj, ...
    'BusyMode', 'drop', ...
	'ExecutionMode', 'fixedSpacing', ...
	'Period', cfg_timer_interval);

% Get the first position
if isempty(pos_ranges)
    cur_pos = [];
else
    cur_pos = pos_ranges(1,:);
end

% Show the first slice
udata = get(fig, 'UserData');
[first_im, first_time] = callback_list.getImage(...
    mov, cur_pos, udata.colormap_spec);

mov_im = image(first_im, 'Parent', mov_ax);
set(mov_ax, 'XTick', []);
set(mov_ax, 'YTick', []);
set(mov_ax, 'Box', 'on');

if ~isempty(first_im)
    set(mov_ax, ...
        'XLim', [1, size(first_im,2)], ...
        'YLim', [1, size(first_im,1)], ...
        'YDir', 'reverse');
else
    set(mov_ax, ...
        'XLim', [0, 100], ...
        'YLim', [0, 100], ...
        'YDir', 'reverse');
end

if ~isempty(first_time)
    set(time_label, 'String', sprintf('Time: %0.3f', first_time));
else
    set(time_label, 'String', '');
end

% Push variables into figure user data
udata.mov = mov;
udata.ch_ids = ch_ids;
udata.callback_list = callback_list;
udata.pos_ranges = pos_ranges;
udata.poly_ctrl = poly_ctrl;
udata.ellipse_ctrl = ellipse_ctrl;
udata.plot_ctrl = plot_ctrl;
udata.play_btn_play_image = play_btn_play_image;
udata.play_btn_pause_image = play_btn_pause_image;
udata.pos_text_labels = pos_text_labels;
udata.pos_play_btns = pos_play_btns;
udata.pos_slider_ctrls = pos_slider_ctrls;
udata.pos_num_edits = pos_num_edits;
udata.rate_label = rate_label;
udata.rate_ctrl = rate_ctrl;
udata.backup_rate = cfg_default_rate;
udata.time_label = time_label;
udata.mov_ax = mov_ax;
udata.mov_im = mov_im;
udata.cmap_fig = cmap_fig;

udata.addListenerToROI = @(uuid,h) addListenerToROI(fig,uuid,h);
udata.removeListenerFromROI = @(uuid,h) removeListenerFromROI(fig,uuid,h);

if ~isempty(first_im)
    udata.aspect_ratio = size(first_im,2) / size(first_im,1);
else
    udata.aspect_ratio = 1;
end

udata.play_pos_idx = [];
udata.play_start_time = [];
udata.play_start_val = [];
udata.timer_obj = timer_obj;
udata.roi_list = [];

set(fig, 'UserData', udata);

% Update the control positions
updatePositions(fig);

% Set the callbacks
for pos_idx=1:size(pos_ranges,2)
    addlistener(pos_slider_ctrls(pos_idx), 'ContinuousValueChange', ...
        @(src,evt) sliderCallback(fig, src, evt));
    
    set(pos_num_edits(pos_idx), 'Callback', ...
        @(src, evt) posEditCallback(fig, src, evt));
    
    set(pos_play_btns(pos_idx), 'Callback', ...
        @(src,evt) posPlayCallback(fig, pos_idx, src, evt));
end

set(rate_ctrl, 'Callback', ...
    @(src, evt) rateCtrlCallback(fig, src, evt));

set(timer_obj, ...
    'TimerFcn', @(src, evt) timerFireCallback(fig, src, evt), ...
    'StartFcn', @(src, evt) timerStartCallback(fig, src, evt), ...
    'StopFcn', @(src, evt) timerStopCallback(fig, src, evt));

set(fig, ...
    'KeyPressFcn', @(src, evt) figKeyPressCallback(fig, src, evt), ...
    'ResizeFcn', @(src, evt) figResizeCallback(fig, src, evt), ...
    'DeleteFcn', @(src, evt) figDeleteCallback(fig, src, evt));
 
set(poly_ctrl, ...
    'ClickedCallback', @(src, evt) polyToolCallback(fig, src, evt));
set(ellipse_ctrl, ...
    'ClickedCallback', @(src, evt) ellipseToolCallback(fig, src, evt));
set(plot_ctrl, ...
    'ClickedCallback', @(src, evt) plotToolCallback(fig, src, evt));
set(cmap_ctrl, ...
    'ClickedCallback', @(src, evt) cmapToolCallback(fig, src, evt));

end

function updatePositions(fig)
% This function recalculates the positions for each of the
% controls based on the current size of the figure.

% Note that the slider height must be 20 because slider height is
% constrained on mac, so this is the only value that give correct
% results

cfg_plabel_hmargin = 15;
cfg_plabel_vmargin = 5;
cfg_plabel_width = 80;

cfg_play_hmargin = 5;
cfg_play_vmargin = 5;
cfg_play_width = 20;
cfg_play_height = 20;

cfg_slider_hmargin = 5;
cfg_slider_vmargin = 2;
cfg_slider_height = 20;
cfg_slider_spacing = 5;

cfg_nedit_hmargin = 5;
cfg_nedit_vmargin = 5;
cfg_nedit_width = 50;
cfg_nedit_height = 20;

cfg_rlabel_hmargin = 15;
cfg_rlabel_vmargin = 14;
cfg_rlabel_width = 65;

cfg_redit_hmargin = 5;
cfg_redit_width = 50;
cfg_redit_height = 20;
cfg_redit_top_margin = 15;

cfg_time_hmargin = 10;
cfg_time_vmargin = 14;
cfg_time_width = 120;

cfg_mov_margins = 20;

% The above configuration was tailored to Mac, but matlab doesn't
% provide proper layout so we have to tweak it so that it doesn't
% look ridiculous on other platforms.  This pretty damned pathetic
% but that is life with matlab.

if isunix
    cfg_slider_vmargin = 5;
    cfg_rlabel_width = 80;
end


fig_pos = get(fig, 'Position');
fig_width = fig_pos(3);
fig_height = fig_pos(4);

udata = get(fig, 'UserData');

% Calculate the font height. MATLAB R2014a does not provide proper
% font metrics so we just have to take a wild guess at the total
% height.
font_height = get(udata.rate_label, 'FontSize') * 1.45;

% Configure the framerate editor
rlabel_xpos = cfg_rlabel_hmargin;
rlabel_ypos = cfg_rlabel_vmargin;
rlabel_width = cfg_rlabel_width;
rlabel_height = font_height;

redit_xpos = rlabel_xpos + rlabel_width + cfg_redit_hmargin;
redit_ypos = rlabel_ypos;
redit_width = cfg_redit_width;
redit_height = cfg_redit_height;

set(udata.rate_label, 'Position', ...
    [rlabel_xpos, rlabel_ypos, rlabel_width, rlabel_height]);
set(udata.rate_ctrl, 'Position', ...
    [redit_xpos, redit_ypos, redit_width, redit_height]);

% Configure the time label
time_xpos = fig_width - cfg_time_width - cfg_time_hmargin;
time_ypos = cfg_time_vmargin;
time_width = cfg_time_width;
time_height = font_height;

set(udata.time_label, 'Position', ...
    [time_xpos, time_ypos, time_width, time_height]);

% Configure the position sliders
cur_ypos = redit_ypos + redit_height + cfg_redit_top_margin;

for pos_idx=size(udata.pos_ranges, 2):-1:1
    % Position label position
    plabel_xpos = cfg_plabel_hmargin;
    plabel_ypos = cfg_plabel_vmargin + cur_ypos;
    plabel_width = cfg_plabel_width;
    plabel_height = font_height;
    
    % Play button position
    play_xpos = plabel_xpos + plabel_width + cfg_play_hmargin;
    play_ypos = cfg_play_vmargin + cur_ypos;
    play_width = cfg_play_width;
    play_height = cfg_play_height;
    
    % Frame edit position
    nedit_xpos = fig_width - cfg_nedit_hmargin - cfg_nedit_width;
    nedit_ypos = cfg_nedit_vmargin + cur_ypos;
    nedit_width = cfg_nedit_width;
    nedit_height = cfg_nedit_height;
    
    % Slider position
    slider_xpos = play_xpos + play_width + cfg_slider_hmargin;
    slider_ypos = cfg_slider_vmargin + cur_ypos;
    slider_width = nedit_xpos - slider_xpos - cfg_slider_hmargin;
    slider_height = cfg_slider_height;
    slider_width = max([1, slider_width]);
    
    set(udata.pos_text_labels(pos_idx), 'Position', ...
        [plabel_xpos, plabel_ypos, plabel_width, plabel_height]);
    set(udata.pos_play_btns(pos_idx), 'Position', ...
        [play_xpos, play_ypos, play_width, play_height]);
    set(udata.pos_slider_ctrls(pos_idx), 'Position', ...
        [slider_xpos, slider_ypos, slider_width, slider_height]);
    set(udata.pos_num_edits(pos_idx), 'Position', ...
        [nedit_xpos, nedit_ypos, nedit_width, nedit_height]);
    
    cur_ypos = cur_ypos + cfg_slider_height + cfg_slider_spacing;
end

% Main movie bounds
mov_bnds_xpos = cfg_mov_margins;
mov_bnds_ypos = slider_ypos + slider_height + cfg_mov_margins;
mov_bnds_width = fig_width - mov_bnds_xpos - cfg_mov_margins;
mov_bnds_height = fig_height - mov_bnds_ypos - cfg_mov_margins;

mov_bnds_width = max([1, mov_bnds_width]);
mov_bnds_height = max([1, mov_bnds_height]);

mov_bnds_aspect_ratio = mov_bnds_width / mov_bnds_height;

if udata.aspect_ratio > mov_bnds_aspect_ratio
    mov_xpos = mov_bnds_xpos;
    mov_width = mov_bnds_width;
    
    mov_ypos = mov_bnds_ypos + 0.5 * mov_bnds_height ...
        - 0.5 * mov_width / udata.aspect_ratio;
    mov_height = mov_width / udata.aspect_ratio;
else
    mov_ypos = mov_bnds_ypos;
    mov_height = mov_bnds_height;
    
    mov_xpos = mov_bnds_xpos + 0.5 * mov_bnds_width ...
        - 0.5 * mov_height * udata.aspect_ratio;
    mov_width = mov_height * udata.aspect_ratio;
end

set(udata.mov_ax, 'Position', [mov_xpos, mov_ypos, ...
    mov_width, mov_height]);
end

function colormapChanged(fig, colormap_spec)
% This function is invoked whenever the colormap is altered

udata = get(fig, 'UserData');

udata.colormap_spec = colormap_spec;

set(fig, 'UserData', udata);

% Find the current position value
pos_val = zeros(1, size(udata.pos_ranges, 2));

for pos_idx=1:size(udata.pos_ranges, 2)
    val = get(udata.pos_slider_ctrls(pos_idx), 'Value');
    pos_val(pos_idx) = round(val);
end

im = udata.callback_list.getImage(udata.mov, ...
    pos_val, udata.colormap_spec);

set(udata.mov_im, 'CData', im);

if ~isempty(im)
    set(udata.mov_ax, ...
        'XLim', [1, size(im,2)], ...
        'YLim', [1, size(im,1)], ...
        'YDir', 'reverse');
else
    set(udata.mov_ax, ...
        'XLim', [0, 100], ...
        'YLim', [0, 100], ...
        'YDir', 'reverse');
end

end

function setPosition(fig, pos)
% This function moves the player to a new position

udata = get(fig, 'UserData');

% Retrive the image
[im,time] = udata.callback_list.getImage(...
    udata.mov, pos, udata.colormap_spec);

if ~isempty(im)
    aspect_ratio = size(im, 2) / size(im, 1);
else
    aspect_ratio = 1;
end

if abs(aspect_ratio - udata.aspect_ratio) > 10*eps
    udata.aspect_ratio = aspect_ratio;
    
    set(fig, 'UserData', udata);
    updatePositions(fig);
end

% Update the time label
if ~isempty(time)
    set(udata.time_label, 'String', sprintf('Time: %0.3f', time));
else
    set(udata.time_label, 'String', '');
end

% Update the image data
set(udata.mov_im, 'CData', im);

if ~isempty(im)
    set(udata.mov_ax, ...
        'XLim', [1, size(im,2)], ...
        'YLim', [1, size(im,1)], ...
        'YDir', 'reverse');
else
    set(udata.mov_ax, ...
        'XLim', [0, 100], ...
        'YLim', [0, 100], ...
        'YDir', 'reverse');
end

set(fig, 'UserData', udata);

end

function pos_val = getCurrentPosition(fig)
% This function retrieves the current position from the sliders

udata = get(fig, 'UserData');

pos_val = zeros(1, size(udata.pos_ranges, 2));
for pos_idx=1:size(udata.pos_ranges, 2)
    pos_val(pos_idx) = round(get(udata.pos_slider_ctrls(pos_idx), 'Value'));
end

end

function sliderCallback(fig, ~, ~)
% This function is invoked when the slider value is changed.

udata = get(fig, 'UserData');

% Stop playback if running
stop(udata.timer_obj);
set(fig, 'UserData', udata);

% Find the current position value
pos_val = zeros(1, size(udata.pos_ranges, 2));

for pos_idx=1:size(udata.pos_ranges, 2)
    min_val = get(udata.pos_slider_ctrls(pos_idx), 'Min');
    max_val = get(udata.pos_slider_ctrls(pos_idx), 'Max');
    val = get(udata.pos_slider_ctrls(pos_idx), 'Value');
    val = min([max([round(val), min_val]), max_val]);
    
    set(udata.pos_num_edits(pos_idx), 'String', num2str(val));
    
    pos_val(pos_idx) = val;
end

% Display the new state
setPosition(fig, pos_val);
end

function posEditCallback(fig, ~, ~)
% This function is invoked when the position text edit changes

udata = get(fig, 'UserData');

% Stop playback if running
stop(udata.timer_obj);
set(fig, 'UserData', udata);

% Find the current position vector
pos_val = zeros(1, size(udata.pos_ranges, 2));
is_changed = false(1, size(udata.pos_ranges, 2));

for pos_idx=1:size(udata.pos_ranges, 2)
    min_val = get(udata.pos_slider_ctrls(pos_idx), 'Min');
    max_val = get(udata.pos_slider_ctrls(pos_idx), 'Max');
    edit_str = get(udata.pos_num_edits(pos_idx), 'String');
    slider_num = get(udata.pos_slider_ctrls(pos_idx), 'Value');
    slider_num = round(slider_num);

    edit_num = str2double(edit_str);
    if isnan(edit_num) || edit_num < min_val || edit_num > max_val
        set(udata.pos_num_edits(pos_idx), 'String', num2str(slider_num));
        edit_num = slider_num;
    else
        edit_num = min([max([round(edit_num), min_val]), max_val]);
        
        if abs(edit_num - slider_num) > 10*eps
            is_changed(pos_idx) = true;
        end
        
        set(udata.pos_slider_ctrls(pos_idx), 'Value', edit_num);
    end

    pos_val(pos_idx) = edit_num;
end

% Give focus to first modified object
first_idx = find(is_changed, 1, 'first');
if ~isempty(first_idx)
    uicontrol(udata.pos_slider_ctrls(first_idx));
end

% Display the current state
setPosition(fig, pos_val);
end

function posPlayCallback(fig, gen_pos_idx, ~, ~)
% This function is invoked when any play button is pushed

udata = get(fig, 'UserData');

running_prop = get(udata.timer_obj, 'Running');
if strcmp(running_prop, 'off')
    % Store the position we should play
    udata.play_pos_idx = gen_pos_idx;
    set(fig, 'UserData', udata);
    
    % Start playback timer
    start(udata.timer_obj);
else
    % We are already running, so do a pause
    stop(udata.timer_obj);
end

end 

function rateCtrlCallback(fig, ~, ~)
% This function is invoked when the framerate edit is modified.

udata = get(fig, 'UserData');

% If the timer is running, then we need to reset the
% start position so that we don't skip frames
running_prop = get(udata.timer_obj, 'Running');
if strcmp(running_prop, 'on')
    udata.play_start_time = clock;
    
    slider = udata.pos_slider_ctrls(udata.play_pos_idx);
    start_val = get(slider, 'Value');
    start_val = round(start_val);
    
    udata.play_start_val = start_val;
    set(fig, 'UserData', udata);
end

val = get(udata.rate_ctrl, 'String');
num = str2double(val);

if isnan(num)
    num = udata.backup_rate;
    set(udata.rate_ctrl, 'String', num2str(num));
else
    udata.backup_rate = num;
    set(fig, 'UserData', udata);
end

% Give focus to the first position slider
uicontrol(udata.pos_slider_ctrls(1));

end

function timerFireCallback(fig, ~, ~)
% This function is invoked with the play timer fires.

udata = get(fig, 'UserData');

slider = udata.pos_slider_ctrls(udata.play_pos_idx);
nedit = udata.pos_num_edits(udata.play_pos_idx);

% Find time elapsed since playback started
if isempty(udata.play_start_time)
    udata.play_start_time = clock;
    
    start_val = get(slider, 'Value');
    start_val = round(start_val);
    
    udata.play_start_val = start_val;
    set(fig, 'UserData', udata);
    
    elapsed = 0.0;
else  
    elapsed = etime(clock, udata.play_start_time);
end

% Find minimum and maximum position values for played position dimension
min_val = get(slider, 'Min');
max_val = get(slider, 'Max');

% Retrieve the current framerate
framerate = get(udata.rate_ctrl, 'String');
framerate = str2double(framerate);
assert(~isnan(framerate), 'Unable to get a valid framerate');

% Calculate the current position with looping
abs_pos = udata.play_start_val + elapsed * framerate;
rel_pos = abs_pos - floor((abs_pos - min_val) / ...
    (max_val - min_val))*(max_val - min_val);

new_pos = min([max([round(rel_pos), min_val]), max_val]);

% Push new position to the appropriate slider and text edit
set(slider, 'Value', new_pos);
set(nedit, 'String', num2str(new_pos));

% Find the new position vector
assem_pos = zeros(1, size(udata.pos_ranges, 2));
for pos_idx=1:size(udata.pos_ranges, 2)
    slider_num = get(udata.pos_slider_ctrls(pos_idx), 'Value');
    assem_pos(pos_idx) = round(slider_num);
end

% Display the new state
setPosition(fig, assem_pos);
end

function timerStartCallback(fig, ~, ~)
% This function is invoked when the play timer starts.

udata = get(fig, 'UserData');

for pos_idx=1:size(udata.pos_ranges, 2)
    if pos_idx ~= udata.play_pos_idx
        set(udata.pos_play_btns(pos_idx), 'Enable', 'off');
    else
        set(udata.pos_play_btns(pos_idx), 'CData', udata.play_btn_pause_image);
    end
    
    set(udata.pos_num_edits(pos_idx), 'Enable', 'off');
end

udata.play_start_time = [];
udata.play_start_pos = [];
set(fig, 'UserData', udata);

end

function timerStopCallback(fig, ~, ~)
% This function is invoked when the play timer stops.

udata = get(fig, 'UserData');
fig_bkg = get(fig, 'Color');

for pos_idx=1:size(udata.pos_ranges, 2)
    set(udata.pos_play_btns(pos_idx), 'Enable', 'on');
    set(udata.pos_play_btns(pos_idx), 'CData', udata.play_btn_play_image);
    set(udata.pos_num_edits(pos_idx), 'Enable', 'on');
    
    % There is a bug in MATLAB R2014a that causes it to render edits that
    % have been renabled with the wrong color, worse still, if we manually
    % overwrite the color it doesn't stick unless it is different then the
    % current color, so we have to set the color twice: once to a random
    % number and then once to the correct number
    set(udata.pos_num_edits(pos_idx), 'BackgroundColor', rand(1, 3));
    set(udata.pos_num_edits(pos_idx), 'BackgroundColor', fig_bkg);
end

end

function figKeyPressCallback(fig, ~, ~)
% This function is invoked when a key is pressed.

udata = get(fig, 'UserData');

% There are currently no recognized key press commands

set(fig,'UserData', udata);
end

function figResizeCallback(fig, ~, ~)
% This function is invoked when the figure is resized.

updatePositions(fig);
end

function figDeleteCallback(fig, ~, ~)
% This function is invoked when the figure is deleted.

udata = get(fig, 'UserData');

delete(udata.cmap_fig);

stop(udata.timer_obj);
delete(udata.timer_obj);
end

function list = addHandleToListeners(h, list)
% This function adds the passed handle object to the cell array list
% if it is not already there

for idx=1:length(list)
    if h == list{idx}
        return;
    end
end

list{end+1} = h;

end

function list = removeHandleFromListeners(h, list)
% This function removes the passed plot object from the cell array
% list.

for idx=1:length(list)
    if h == list{idx}
        list(idx) = [];
    end
end

end

function addListenerToROI(fig, uuid, h)
% This function adds an ROI to the listener list. 

udata = get(fig, 'UserData');
roi_list = udata.roi_list;

for roi_idx=1:length(roi_list)
    if strcmp(uuid, roi_list(roi_idx).uuid)
        roi_list(roi_idx).listeners = addHandleToListeners(...
            h, roi_list(roi_idx).listeners);
        
        udata.roi_list = roi_list;
        set(fig, 'UserData', udata);
        handleROIChange(fig, roi_list(roi_idx).handle);
    end
end

end

function removeListenerFromROI(fig, uuid, h)
% This function removes an ROI from the listener list.

udata = get(fig, 'UserData');
roi_list = udata.roi_list;

for roi_idx=1:length(roi_list)
    if strcmp(uuid, roi_list(roi_idx).uuid)
        roi_list(roi_idx).listeners = removeHandleFromListeners(...
            h, roi_list(roi_idx).listeners);
    end
end

udata.roi_list = roi_list;
set(fig, 'UserData', udata);

end

function replaceROIContextMenu(fig, h)
% This function replaces the context menu current associated
% with the passed ROI handle 'h'

cfg_color_choices = { ...
    {'Blue', [0, 0, 1]}, ...
    {'Red', [1, 0, 0]}, ...
    {'Green' [0, 1, 0]}, ...
    {'Yellow', [1, 1, 0]}, ...
    {'Magenta', [1, 0, 1]}, ...
    {'Cyan', [0, 1, 1]}, ...
    {'Black', [0, 0, 0]}...
    };

% Create the main context menu
cmenu = uicontextmenu('Parent', fig);

% Create context menu for setting color
color_menu = uimenu(cmenu, ...
	'Label', 'Set Color');

color_submenus = zeros(1, numel(cfg_color_choices));
for submenu_idx=1:numel(cfg_color_choices)
	color_submenus(submenu_idx) = uimenu(color_menu,...
		'Label', cfg_color_choices{submenu_idx}{1}, ...
		'Callback', @(src, evt) roiColorCallback(fig, h, ...
        cfg_color_choices{submenu_idx}{2}, src, evt));
end

% Create context menu for deleting the roi
uimenu(cmenu, ...
    'Label', 'Delete', ...
    'Callback', @(src,evt) roiDeleteCallback(fig, h, src, evt));

% Assign the context menu to the object
cmenu_obj = findobj(h, 'Type', 'line', '-or', 'Type', 'patch');  
	set(cmenu_obj, 'uicontextmenu', cmenu);
       
set(cmenu_obj, 'uicontextmenu', cmenu);       

end

function handleROIChange(fig, roi)
% This function updates the passed ROI in response to a change
% in shape or position

udata = get(fig, 'UserData');

for idx=1:length(udata.roi_list)
    if udata.roi_list(idx).handle == roi
        
        listeners = udata.roi_list(idx).listeners;
        if length(listeners) >= 1
            trimmed_roi_list = rmfield(udata.roi_list(idx), ...
                {'uuid', 'listeners'});
            
            dset = udata.callback_list.processROIs(...
                udata.mov, trimmed_roi_list);
        end
        
        listener_idx = 1;
        while listener_idx <= length(listeners)
            if ~ishghandle(listeners{listener_idx})
                % Remove hanging listener handle     
                listeners(listener_idx) = [];
            else
                % Push new data to lineseries object
                set(listeners{listener_idx}, ...
                    'XData', dset{1}(1,:), ...
                    'YData', dset{1}(2,:));
            end
            
            listener_idx = listener_idx + 1;
        end
        
        udata.roi_list(idx).listeners = listeners;
        
        break;
    end
end

set(fig, 'UserData', udata);

end

function addROIToFig(fig, roi)
% This function adds the roi to the list of ROIs that
% are kept by the figure;

udata = get(fig, 'UserData');

addNewPositionCallback(roi, ...
    @(pos) roiChangedCallback(fig, roi, pos));

udata.roi_list(end+1).uuid = char(java.util.UUID.randomUUID);
udata.roi_list(end).handle = roi;
udata.roi_list(end).channel = udata.ch_ids(1);
udata.roi_list(end).home_pos = getCurrentPosition(fig);
udata.roi_list(end).scan_pos = 1;
udata.roi_list(end).listeners = {};

set(fig, 'UserData', udata);
end


function removeROIFromFig(fig, roi)
% This removes the passed ROI from the figure.

udata = get(fig, 'UserData');

match_idx = find([udata.roi_list.handle] == roi, 1, 'first');
if isempty(match_idx)
    return;
end

old_roi = udata.roi_list(match_idx).handle;
udata.roi_list(match_idx) = [];

delete(old_roi);

if ~isempty(udata.roi_list)
    set(udata.plot_ctrl, 'Enable', 'on');
else
    set(udata.plot_ctrl, 'Enable', 'off');
end

set(fig,'UserData', udata);
end


function roiColorCallback(~, roi, color, ~, ~)

roi.setColor(color);
% TODO: update the linked figures

end

function roiChangedCallback(fig, roi, ~)
% This function is invoked whenever the roi is moved or resized.

handleROIChange(fig, roi);
end

function roiDeleteCallback(fig, roi, ~, ~)
% This function is invoked whenever the roi should be deleted.

removeROIFromFig(fig, roi);
end


function disableToolbar(fig)
% This function disables the toolbar buttons

udata = get(fig, 'UserData');
set(udata.poly_ctrl, 'Enable', 'off');
set(udata.ellipse_ctrl, 'Enable', 'off');
set(udata.plot_ctrl, 'Enable', 'off');
end

function enableToolbar(fig)
% This function enables the toolbar buttons

udata = get(fig, 'UserData');
set(udata.poly_ctrl, 'Enable', 'on');
set(udata.ellipse_ctrl, 'Enable', 'on');

if ~isempty(udata.roi_list)
    set(udata.plot_ctrl, 'Enable', 'on');
else
    set(udata.plot_ctrl, 'Enable', 'off');
end

end

function polyToolCallback(fig, ~, ~)
% This function is invoked when the polygon tool is invoked.

udata = get(fig, 'UserData');
disableToolbar(fig);

h = impoly(udata.mov_ax);
replaceROIContextMenu(fig, h);

if ~isempty(h)
    fcn = makeConstrainToRectFcn('impoly', ...
        get(udata.mov_ax,'XLim'), get(udata.mov_ax,'YLim'));
    setPositionConstraintFcn(h, fcn);

    addROIToFig(fig, h);
end

enableToolbar(fig);

end

function ellipseToolCallback(fig, ~, ~)
% This function is invoked when the ellipse tool is invoked.

udata = get(fig, 'UserData');
disableToolbar(fig);

h = imellipse(udata.mov_ax);
replaceROIContextMenu(fig, h);

if ~isempty(h)
    fcn = makeConstrainToRectFcn('imellipse', ...
        get(udata.mov_ax,'XLim'), get(udata.mov_ax,'YLim'));
    setPositionConstraintFcn(h, fcn);

    addROIToFig(fig, h);
end

enableToolbar(fig);

end

function plotToolCallback(fig, ~, ~)
% This function is invoked when the plot tool is invoked.

cfg_trace_fig_xpos = 50;
cfg_trace_fig_ypos = -50;
cfg_trace_fig_width = 400;
cfg_trace_fig_height = 150;

udata = get(fig, 'UserData');
fig_pos = get(fig, 'Position');

trimmed_roi_list = rmfield(udata.roi_list, ...
	{'uuid', 'listeners'});
dset = udata.callback_list.processROIs(udata.mov, trimmed_roi_list);

colors = cell(1, length(udata.roi_list));
live_info.creator = fig;
live_info.uuids = cell(1, length(udata.roi_list));

for idx=1:length(udata.roi_list)
    colors{idx} = getColor(udata.roi_list(idx).handle);
    live_info.uuids{idx} = udata.roi_list(idx).uuid;
end

trace_fig_height = cfg_trace_fig_height * length(dset);
trace_fig_height = max([256, trace_fig_height]);
trace_fig_height = min([1024, trace_fig_height]);

simul_fig = figure('Visible', 'off');
set(simul_fig, 'Position', [fig_pos(1) + cfg_trace_fig_xpos, ...
    fig_pos(2) + cfg_trace_fig_ypos, cfg_trace_fig_width, ...
    trace_fig_height]);
set(simul_fig, 'Visible', 'on');

nia_displaySimulPlots(simul_fig, dset, 'Position', 'Intensity', ...
    colors, live_info);

end

function cmapToolCallback(fig, ~, ~)
% This function is invoked when the cmap button on the 
% toolbar is clicked.

udata = get(fig, 'UserData');

if ishghandle(udata.cmap_fig)
    figure(udata.cmap_fig);
end

end

