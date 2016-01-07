function colormap_spec = nia_selectColormap(fig, ch_ids, hist_info, callback)
%NIA_SELECTCOLORMAP Open GUI to select colormaps
%   colormap_spec = nia_select_colormap(fig, hist_info, callback) opens a
%   graphical window that displays the histogram for a given image, and
%   provides the user with graphical controls that can be used to select a
%   colormap and manipulate the colormap range. It returns a colormap spec
%   based on the input values. If a callback function is provided then it
%   is assumed that this function is managed, as such it will hide when it
%   receives a close request rather delete itself. The managing figure must
%   delete it when appropriate.
%
%   This function accepts the following arguments:
%
%       fig - Figure to use for axes.  If empty a new figure
%           will be created
%
%       ch_ids - A [1xN] array containing the channel identifier
%           for each channel in hist_info
%
%       hist_info - A structure containing information about
%           the image histogram, more information can be found
%           below.
%
%       callback - An optional function handle that is invoked
%           and passed a colormap specification whenever the
%           current colormap is altered.  The colormap
%           specification is described below.
%
%   The argument hist_info specifies the histogram information
%   that informs selection of the colormap.  It must be a structure
%   array with an element for each channel and the fields:
%
%       min - Minimum intensity value
%       max - Maximum intensity value
%       dist - A 2xN array specifying the histogram counts. The
%           first row contains the count location, and the 
%           second row contains the count value. This
%           field may also be left empty, in which case the
%           histogram plot is left empty
%
%   The colormap specification passed to the callback function
%   is a structure array with the fields:
%
%       min - A scalar value indicating minimum value for colormap
%       max - A scalar value indicating maximum value for colormap
%       colormap - A Nx3 matrix for use a colormap
%       nan_color - A 1x3 matrix used for NaN
%

cfg_label_font_size = 12;

% Fill in optional arguments
if nargin < 2
    callback = [];
end

% Check input arguments
if ~ishghandle(fig)
    error 'The argument ''fig'' must be a figure handle';
end

if ~isfloat(ch_ids) || ~isreal(ch_ids) || ~isrow(ch_ids)
    error 'The argument ''ch_ids'' has an invalid type';
end

if nnz(abs(ch_ids - round(ch_ids))) > 0
    error 'The argument ''ch_ids'' must have integer values';
end

if ~isstruct(hist_info) || ~isrow(hist_info)
    error 'The argument ''hist_info'' has an invalid type';
end

hist_info_fields = {'min', 'max', 'dist'};
[hist_info_ok, hist_info_msg] = nia_hasValidFieldNames(...
    hist_info, hist_info_fields, hist_info_fields);
if ~hist_info_ok
    error(hist_info_msg, 'hist_info');
end

for idx=1:length(hist_info)
    if ~isfloat(hist_info(idx).min) || ~isreal(hist_info(idx).min) || ...
            ~isscalar(hist_info(idx).min)
        error 'The argument ''hist_info.min'' has an invalid type';
    end
    
    if ~isfloat(hist_info(idx).max) || ~isreal(hist_info(idx).max) || ...
            ~isscalar(hist_info(idx).max)
        error 'The argument ''hist_info.max'' has an invalid type';
    end
    
    if ~isfloat(hist_info(idx).dist) || ~isreal(hist_info(idx).dist) || ...
            ~ismatrix(hist_info(idx).dist)
        error 'The argument ''hist_info.dist'' has an invalid type';
    end
    
    if ~isempty(hist_info(idx).dist)
        if size(hist_info(idx).dist, 1) ~= 2
            error 'The argument ''hist_info.dist'' has invalid dimensions';
        end
        
        if nnz(hist_info(idx).dist(2,:) < 0) > 0
            error 'The argument ''hist_info.dist'' has invalid values';
        end
    end
end

if length(hist_info) ~= length(ch_ids)
    error 'The length of ''ch_ids'' and ''hist_info'' must match';
end

if ~isempty(callback)
    if ~isa(callback, 'function_handle')
        error 'The argument ''callback'' has an invalid type';
    end
end

% Create the figure if necessary
if isempty(fig)
    fig = figure;
end

set(fig, 'MenuBar', 'none');
fig_bkg = get(fig, 'Color');

% Create the colormap choices, these must be aligned
% with switch statement in get_colormap_spec()
cmap_choices = {'uni_colormap1', 'uni_colormap2', ...
    'red', 'green', 'blue', 'yellow', 'magenta', ...
    'cyan', 'darkfield', 'brightfield'};
cmap_start_rotate = 3;
cmap_end_rotate = 6;

% Create the channel selector
cur_ch = 1;
chan_num = length(hist_info);
chan_descr_list = arrayfun(@(x) num2str(x), ch_ids, ...
    'UniformOutput', false);

chan_text = uicontrol(...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'String', 'Channel', ...
    'FontSize', cfg_label_font_size, ...
    'Background', fig_bkg);
chan_ctrl = uicontrol(...
    'Style', 'popupmenu', ...
    'Units', 'pixels', ...
    'String', chan_descr_list);

% Set default values, preallocate array
chan_list(chan_num).min = [];
chan_list(chan_num).max = [];
chan_list(chan_num).cmap = [];

for chan_idx=1:chan_num
    chan_list(chan_idx).min = hist_info(chan_idx).min;
    chan_list(chan_idx).max = hist_info(chan_idx).max;
    
    cmap_idx = cmap_start_rotate + ...
        mod(chan_idx-1, cmap_end_rotate-cmap_start_rotate+1);
    chan_list(chan_idx).cmap = cmap_idx;
end

% Create colormap selector
cmap_text = uicontrol(...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'String', 'Colormap', ...
    'FontSize', cfg_label_font_size, ...
    'Background', fig_bkg);
cmap_ctrl = uicontrol(...
    'Style', 'popupmenu', ...
    'Units', 'pixels', ...
    'String', cmap_choices, ...
    'Value', chan_list(cur_ch).cmap);

% Create the minimum and maximum sliders
min_text = uicontrol(...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'String', 'Minimum', ...
    'FontSize', cfg_label_font_size, ...
    'Background', fig_bkg);
min_ctrl = uicontrol(...
    'Style', 'slider', ...
    'Units', 'pixels', ...
    'Value', chan_list(cur_ch).min, ...
    'Min', hist_info(cur_ch).min, ...
    'Max', hist_info(cur_ch).max);

max_text = uicontrol(...
    'Style', 'text', ...
    'Units', 'pixels', ...
    'String', 'Maximum', ...
    'FontSize', cfg_label_font_size, ...
    'Background', fig_bkg);
max_ctrl = uicontrol(...
    'Style', 'slider', ...
    'Units', 'pixels', ...
    'Value', chan_list(cur_ch).max, ...
    'Min', hist_info(cur_ch).min, ...
    'Max', hist_info(cur_ch).max);

% Create the histogram axes
hist_ax = axes;

if ~isempty(hist_info(cur_ch).dist)
    bar_h = bar(hist_ax, hist_info(cur_ch).dist(1,:), ...
        hist_info(cur_ch).dist(2,:), 1.0);
    
    set(bar_h, 'FaceColor', [0.5, 0.5, 0.5]);
else
    bar_h = [];
end

set(hist_ax, ...
    'Unit', 'pixels', ...
    'Box', 'on', ...
    'XLim', [hist_info(cur_ch).min, ...
            hist_info(cur_ch).max]);
xlabel(hist_ax, 'Intensity');
ylabel(hist_ax, 'Count');

% Force Y bottom to zero
hist_ylims = zeros(1, 2);
hist_ylims(1) = 0;
hist_ylims(2) = getPlotUpperBnd(hist_info(cur_ch).dist(2,:));
set(hist_ax, 'YLim', hist_ylims);

% Draw minimum and maximum lines
min_line = line(...
    [chan_list(cur_ch).min, chan_list(cur_ch).min], ...
    [hist_ylims(1), hist_ylims(2)], ...
    'Parent', hist_ax);

set(min_line, 'Color', 'black', 'LineWidth', 2);

max_line = line(...
    [chan_list(cur_ch).max, chan_list(cur_ch).max], ...
    [hist_ylims(1), hist_ylims(2)], ...
    'Parent', hist_ax);

set(max_line, 'Color', 'black', 'LineWidth', 2);

min_max_line = line(...
    [chan_list(cur_ch).min, chan_list(cur_ch).max], ...
    [hist_ylims(1), hist_ylims(2)], ...
    'Parent', hist_ax);

set(min_max_line, 'Color', 'black', 'LineWidth', 2);


% Push data in user data
udata.hist_info = hist_info;
udata.callback = callback;
udata.chan_list = chan_list;
udata.chan_text = chan_text;
udata.chan_ctrl = chan_ctrl;
udata.cmap_text = cmap_text;
udata.cmap_ctrl = cmap_ctrl;
udata.min_text = min_text;
udata.min_ctrl = min_ctrl;
udata.min_line = min_line;
udata.max_text = max_text;
udata.max_ctrl = max_ctrl;
udata.max_line = max_line;
udata.min_max_line = min_max_line;
udata.hist_ax = hist_ax;
udata.bar_h = bar_h;

updatePositions(fig, udata);

set(fig, 'UserData', udata);

% Add callbacks
set(chan_ctrl, 'Callback', @(src, evt) channelCallback(fig, src, evt));
set(cmap_ctrl, 'Callback', @(src, evt) cmapCallback(fig, src, evt));

addlistener(min_ctrl, 'ContinuousValueChange', ...
    @(src,evt) minSliderCallback(fig, src, evt));
addlistener(max_ctrl, 'ContinuousValueChange', ...
    @(src,evt) maxSliderCallback(fig, src, evt));

set(fig, ...
    'ResizeFcn', @(src, evt) figResizeCallback(fig, src, evt), ...
    'CloseRequestFcn', @(src, evt) figCloseRequestCallback(fig, src, evt));

colormap_spec = retrieveColormapSpec(fig);

end


function colormap_spec = retrieveColormapSpec(fig)
% This function retrieves the current colormap spec.

cfg_hist_nelem = 100;

udata = get(fig, 'UserData');

chan_list = udata.chan_list;
num_ch = length(chan_list);

% preallocate colormap_spec array
colormap_spec(num_ch).min = [];
colormap_spec(num_ch).max = [];
colormap_spec(num_ch).colormap = [];
colormap_spec(num_ch).nan_color = [];

for ch_idx=1:num_ch
    colormap_spec(ch_idx).min = chan_list(ch_idx).min;
    colormap_spec(ch_idx).max = chan_list(ch_idx).max;
    
    % These must be aligned with cmap_choices in main
    % function.
    switch chan_list(ch_idx).cmap
        case 1
            % uni1 colormap
            cmap = nia_colormapUni1;
            nan_color = [1, 1, 1];
        case 2
            % uni2 colormap
            cmap = nia_colormapUni2;
            nan_color = [1, 1, 1];
        case 3
            % red colormap
            cmap = [linspace(0, 1, cfg_hist_nelem)', ...
                    zeros(cfg_hist_nelem, 1), ...
                    zeros(cfg_hist_nelem, 1)];
            nan_color = [1, 1, 1];
        case 4
            % green colormap
            cmap = [zeros(cfg_hist_nelem, 1), ...
                    linspace(0, 1, cfg_hist_nelem)', ...
                    zeros(cfg_hist_nelem, 1)];
            nan_color = [1, 1, 1];
        case 5
            % blue colormap
            cmap = [zeros(cfg_hist_nelem, 1), ...
                    zeros(cfg_hist_nelem, 1), ...
                    linspace(0, 1, cfg_hist_nelem)'];
            nan_color = [1, 1, 1];
        case 6
             % yellow colormap
            cmap = [linspace(0, 1, cfg_hist_nelem)', ...
                    linspace(0, 1, cfg_hist_nelem)', ...
                    zeros(cfg_hist_nelem, 1)];
            nan_color = [1, 1, 1];
        case 7
            % magenta colormap
            cmap = [linspace(0, 1, cfg_hist_nelem)', ...
                    zeros(cfg_hist_nelem, 1), ...
                    linspace(0, 1, cfg_hist_nelem)'];
            nan_color = [1, 1, 1];
        case 8
            % cyan colormap
            cmap = [zeros(cfg_hist_nelem, 1), ...
                    linspace(0, 1, cfg_hist_nelem)', ...
                    linspace(0, 1, cfg_hist_nelem)'];
            nan_color = [1, 1, 1];            
        case 9
            % darkfield colormap
            cmap = repmat(linspace(0, 1, cfg_hist_nelem)', 1, 3);
            nan_color = [1, 1, 0];
        case 10
            % brightfield colormap
            cmap = repmat(linspace(1, 0, cfg_hist_nelem)', 1, 3);
            nan_color = [1, 1, 0];   
    end
    
    colormap_spec(ch_idx).colormap = cmap;
    colormap_spec(ch_idx).nan_color = nan_color;
end

end

function updatePositions(fig, udata)
% This function updates the positions of the controls
% based on the current figure dimensions.

cfg_text_right_spacing = 10;
cfg_text_width = 80;

% The slider height must be 20 since this is the only value
% supported (correctly) on Mac as of R2014a
cfg_slider_height = 20;

cfg_slider_voffset = -2;
cfg_chan_width = 150;
cfg_chan_height = 20;
cfg_cmap_width = 150;
cfg_cmap_height = 20;
cfg_ctrl_margins = 10;
cfg_ctrl_height = 20;
cfg_ctrl_spacing = 8;
cfg_hist_margins = 15;
cfg_hist_spacing = 15;

fig_props = get(fig);
fig_pos = fig_props.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);

% Calculate the font height. MATLAB R2014a does not provide proper
% font metrics so we just have to take a wild guess at the baseline
% position.
font_height = get(udata.cmap_text, 'FontSize') * 1.45;

avail_height = fig_height - cfg_ctrl_margins - cfg_hist_margins ...
    - cfg_hist_spacing - 3*cfg_ctrl_spacing - 4*cfg_ctrl_height;

% Configure the histogram axes
hist_xpos = cfg_hist_margins;
hist_ypos = cfg_hist_margins;
hist_width = fig_width - 2*cfg_hist_margins;
hist_height = avail_height;

% Adjust for tight inset
hist_inset = get(udata.hist_ax, 'TightInset');

hist_xpos = hist_xpos + hist_inset(1);
hist_ypos = hist_ypos + hist_inset(2);
hist_width = hist_width - hist_inset(1) - hist_inset(3);
hist_height = hist_height - hist_inset(2) - hist_inset(4);

hist_width = max([1, hist_width]);
hist_height = max([1, hist_height]);

set(udata.hist_ax, 'Position', [hist_xpos, ...
    hist_ypos, hist_width, hist_height]);

% Configure the maximum slider
max_text_xpos = hist_xpos;
max_text_ypos = cfg_hist_margins + avail_height + cfg_hist_spacing;
max_text_width = cfg_text_width;
max_text_height = cfg_ctrl_height;

set(udata.max_text, 'Position', [max_text_xpos, max_text_ypos, ...
    max_text_width, max_text_height]);

max_ctrl_xpos = max_text_xpos + max_text_width + cfg_text_right_spacing;
max_ctrl_ypos = max_text_ypos + cfg_ctrl_height - font_height + cfg_slider_voffset;
max_ctrl_width = hist_xpos + hist_width - max_ctrl_xpos - cfg_ctrl_margins;
max_ctrl_height = cfg_slider_height;

max_ctrl_width = max([1, max_ctrl_width]);
max_ctrl_height = max([1, max_ctrl_height]);

set(udata.max_ctrl, 'Position', [max_ctrl_xpos, max_ctrl_ypos, ...
    max_ctrl_width, max_ctrl_height]);

% Configure the minimum slider
min_text_xpos = max_text_xpos;
min_text_ypos = max_text_ypos + max_text_height + cfg_ctrl_spacing;
min_text_width = cfg_text_width;
min_text_height = cfg_ctrl_height;

set(udata.min_text, 'Position', [min_text_xpos, min_text_ypos, ...
    min_text_width, min_text_height]);

min_ctrl_xpos = min_text_xpos + min_text_width + cfg_text_right_spacing;
min_ctrl_ypos = min_text_ypos + cfg_ctrl_height - font_height + cfg_slider_voffset;
min_ctrl_width = max_ctrl_width;
min_ctrl_height = cfg_slider_height;

min_ctrl_width = max([1, min_ctrl_width]);
min_ctrl_height = max([1, min_ctrl_height]);

set(udata.min_ctrl, 'Position', [min_ctrl_xpos, min_ctrl_ypos, ...
    min_ctrl_width, min_ctrl_height]);

% Configure the colormap controls
cmap_text_xpos = max_text_xpos;
cmap_text_ypos = min_text_ypos + min_text_height + cfg_ctrl_spacing;
cmap_text_width = cfg_text_width;
cmap_text_height = cfg_ctrl_height;

set(udata.cmap_text, 'Position', [cmap_text_xpos, cmap_text_ypos, ...
    cmap_text_width, cmap_text_height]);

cmap_ctrl_xpos = cmap_text_xpos + cmap_text_width + cfg_text_right_spacing;
cmap_ctrl_ypos = cmap_text_ypos + cfg_ctrl_height - font_height;
cmap_ctrl_width = cfg_cmap_width;
cmap_ctrl_height = cfg_cmap_height;

set(udata.cmap_ctrl, 'Position', [cmap_ctrl_xpos, cmap_ctrl_ypos, ...
    cmap_ctrl_width, cmap_ctrl_height]);

% Configure channel controls
chan_text_xpos = max_text_xpos;
chan_text_ypos = cmap_text_ypos + cmap_text_height + cfg_ctrl_spacing;
chan_text_width = cfg_text_width;
chan_text_height = cfg_ctrl_height;

set(udata.chan_text, 'Position', [chan_text_xpos, chan_text_ypos, ...
    chan_text_width, chan_text_height]);

chan_ctrl_xpos = chan_text_xpos + chan_text_width + cfg_text_right_spacing;
chan_ctrl_ypos = chan_text_ypos + cfg_ctrl_height - font_height;
chan_ctrl_width = cfg_chan_width;
chan_ctrl_height = cfg_chan_height;

set(udata.chan_ctrl, 'Position', [chan_ctrl_xpos, chan_ctrl_ypos, ...
    chan_ctrl_width, chan_ctrl_height]);

end

function figResizeCallback(fig, ~, ~)
% This function is invoked when the figure is resized.

udata = get(fig, 'UserData');
updatePositions(fig, udata)
end

function figCloseRequestCallback(fig, ~, ~)
% This function is invoked when the figure receives a close request.

udata = get(fig, 'UserData');

if ~isempty(udata.callback)
    % Hide rather than close the window
    set(fig, 'Visible', 'off');
end

end

function out = getPlotUpperBnd(x)
% This function finds a nice round number for use as an upper
% bound for a plot of the values in x

    x_max = max(x);
    
    if x_max < 0
        error 'Function only implemented for positive numbers';
    elseif x_max == 0
        out = 1;
        return;
    end
    
    good_numbers = [1.0, 1.25, 1.5, 1.75, 2, 3, 4, 5, 6, 7, 8, 9, 10.0];
    
    base = 10^(floor(log10(x_max)));
    trial_numbers = base * good_numbers;
    best_number = trial_numbers(find(trial_numbers > x_max, 1, 'first'));
    
    out = best_number;
end

function channelCallback(fig, ~, ~)
% This function is invoked when the channel selector is altered

udata = get(fig, 'UserData');
new_chan = get(udata.chan_ctrl, 'Value');

set(udata.cmap_ctrl, ...
    'Value', udata.chan_list(new_chan).cmap);

set(udata.min_ctrl, ...
    'Value', udata.chan_list(new_chan).min, ...
    'Min', udata.hist_info(new_chan).min, ...
    'Max', udata.hist_info(new_chan).max);

set(udata.max_ctrl, ...
    'Value', udata.chan_list(new_chan).max, ...
    'Min', udata.hist_info(new_chan).min, ...
    'Max', udata.hist_info(new_chan).max);

hist_ylims = zeros(1, 2);
hist_ylims(1) = 0;
hist_ylims(2) = getPlotUpperBnd(udata.hist_info(new_chan).dist(2,:));
set(udata.hist_ax, 'YLim', hist_ylims);

set(udata.min_line, ...
    'XData', [udata.chan_list(new_chan).min, udata.chan_list(new_chan).min], ...
    'YData', hist_ylims);

set(udata.max_line, ...
    'XData', [udata.chan_list(new_chan).max, udata.chan_list(new_chan).max], ...
    'YData', hist_ylims);

set(udata.min_max_line, ...
    'XData', [udata.chan_list(new_chan).min, udata.chan_list(new_chan).max], ...
    'YData', hist_ylims);

set(udata.hist_ax, 'XLim', [udata.hist_info(new_chan).min, ...
    udata.hist_info(new_chan).max]);

set(udata.bar_h, ...
    'XData', udata.hist_info(new_chan).dist(1,:), ...
    'YData', udata.hist_info(new_chan).dist(2,:));

end

function cmapCallback(fig, ~, ~)
% This function is invoked when the colormap selector is altered.

udata = get(fig, 'UserData');

cur_chan = get(udata.chan_ctrl, 'Value');
cmap_val = get(udata.cmap_ctrl, 'Value');

udata.chan_list(cur_chan).cmap = cmap_val;

set(fig, 'UserData', udata);

if ~isempty(udata.callback)
    udata.callback(retrieveColormapSpec(fig));
end

end

function minSliderCallback(fig, ~, ~)
% This function is invoked when the minimum slider is altered.

udata = get(fig, 'UserData');

min_slider_val = get(udata.min_ctrl, 'Value');
max_slider_val = get(udata.max_ctrl, 'Value');

set(udata.min_line, 'XData', [min_slider_val, min_slider_val]);
set(udata.min_max_line, 'XData', [min_slider_val, max_slider_val]);

cur_chan = get(udata.chan_ctrl, 'Value');
udata.chan_list(cur_chan).min = min_slider_val;

set(fig, 'UserData', udata);

if ~isempty(udata.callback)
    udata.callback(retrieveColormapSpec(fig));
end

end

function maxSliderCallback(fig, ~, ~)
% This function is invoked when the maximum slider is altered.

udata = get(fig, 'UserData');

min_slider_val = get(udata.min_ctrl, 'Value');
max_slider_val = get(udata.max_ctrl, 'Value');

set(udata.max_line, 'XData', [max_slider_val, max_slider_val]);
set(udata.min_max_line, 'XData', [min_slider_val, max_slider_val]);

cur_chan = get(udata.chan_ctrl, 'Value');
udata.chan_list(cur_chan).max = max_slider_val;

set(fig, 'UserData', udata);

if ~isempty(udata.callback)
    udata.callback(retrieveColormapSpec(fig));
end

end
