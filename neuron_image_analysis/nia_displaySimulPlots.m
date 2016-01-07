function nia_displaySimulPlots(fig, data, x_label, y_label, colors, live_info)
%NIA_DISPLAYSIMULPLOTS Show a set of time synched plots
%   nia_displaySimulPlots(fig, data, x_label, y_label, colors, live_info)
%   displays the time series contained in the argument data in the passed
%   figure with the x axes aligned and locked to each other.
%
%   This function accepts the following arguments:
%
%       fig - Figure to use for axes.  If empty a new figure
%           will be created
%
%       data - Data to plot. Data must be in the form of a 
%           1-dimesional cell array, in which every element
%           of the cell array contains the data for a single
%           trace. The data for each trace must be stored as
%           a two-dimensional matrix with two rows, the X
%           values must be stored in the first row, and the
%           Y values must be stored in the second row.
%
%       x_label - Label to use for x axis
%
%       y_label - Label to use for y axes
%
%       colors - Optional argument specifying the color for each
%           trace.  If provided, it must be a cell array with the
%           same dimensions as data, and every element must be
%           valid ColorSpec color. It omitted, all traces are
%           drawn in gray
%
%       live_info - Option argument specifying information necessary
%           to provide live updates.  The argument must be a scalar
%           structure with the fields creator and uuids.  The creator
%           field must be a figure handle.  The uuids field must be
%           a cell array with the same dimensions as data, where 
%           each element is a string that uniquely identifies the
%           data for that trace.
%
%   Example:
%
%       x = linspace(0, 6*pi, 1000);
%       y1 = sin(2*pi*x);
%       y2 = cos(2*pi*x);
%       data = {[x; y1], [x; y2]};
%       colors = {'green', 'blue'};
%       nia_display_simul_plots([], data, 'Time (s)', 'Response', colors);


% Check input arguments

if ~ishghandle(fig)
    error 'The argument ''fig'' has invalid type';
end

if ~iscell(data) || ~isvector(data)
    error 'The argument ''data'' has an invalid type';
end

for idx=1:length(data)
    dset = data{idx};
    
    if ~isfloat(dset) || ~isreal(dset) || ~ismatrix(dset) || ...
            size(dset,1) ~= 2
        error 'The argument ''data'' has invalid elements';
    end
end

if nargin <= 4
    colors = [];   
else
    if ~iscell(colors) || ~isvector(colors) || ...
            length(colors) ~= length(data)
        error 'The argument ''colors'' has an invalid type';
    end
    
    for idx=1:length(colors)
        if ~nia_isColor(colors{idx})
            error 'The argument ''colors'' has an invalid entry';
        end
    end
end

if nargin <= 5
    live_info = [];
else
    fnames = {'creator', 'uuids'};
    
    if ~isstruct(live_info) || length(live_info) ~= 1
        error 'The argument ''live_info'' has an invalid type';
    end
    
    [live_info_ok, live_info_msg] = nia_hasValidFieldNames(...
        live_info, fnames, fnames);
   	if ~live_info_ok
        error(live_info_msg, 'live_info');
    end
    
    if ~ishghandle(live_info.creator)
        error 'The argument ''live_info.creator'' has an invalid type';
    end
    
    if ~iscell(live_info.uuids) || ~isvector(live_info.uuids) || ...
            length(live_info.uuids) ~= length(data)
        error 'The argument ''live_info.uuids'' has an invalid type';
    end
    
    for idx=1:length(live_info.uuids)
        if ~nia_isString(live_info.uuids{idx})
            error 'The argument ''live_info.uuids'' has an invalid entry';
        end
    end
end


% If figure does not exist, then create one.
if isempty(fig)
    fig = figure;
end

udata = [];

% Configure the figure
set(fig, 'MenuBar', 'none');

% Configure the toolbar
% NOTE: handle visibility is set to off below because MATLAB tutorial
% says to do that.  I have no idea why, matlab has a pretty stupid
% event propagation model.
tbar = uitoolbar(fig);

zoom_orig_image = imread('nia_zoom_orig_icon.png');
zoom_orig_ctrl = uipushtool(tbar, ...
    'CData', zoom_orig_image, ...
    'TooltipString', 'Restore original plot dimensions', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

zoom_xy_image = imread('nia_zoom_xy_icon.png');
zoom_xy_ctrl = uitoggletool(tbar, ...
    'CData', zoom_xy_image, ...
    'TooltipString', 'Zoom in both X & Y directions', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

zoom_x_image = imread('nia_zoom_x_icon.png');
zoom_x_ctrl = uitoggletool(tbar, ...
    'CData', zoom_x_image, ...
    'TooltipString', 'Zoom in X direction only', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

zoom_y_image = imread('nia_zoom_y_icon.png');
zoom_y_ctrl = uitoggletool(tbar, ...
    'CData', zoom_y_image, ...
    'TooltipString', 'Zoom in y direction only', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

pan_xy_image = imread('nia_pan_xy_icon.png');
pan_xy_ctrl = uitoggletool(tbar, ...
    'CData', pan_xy_image, ...
    'TooltipString', 'Pan in both X & Y directions', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

pan_x_image = imread('nia_pan_x_icon.png');
pan_x_ctrl = uitoggletool(tbar, ...
    'CData', pan_x_image, ...
    'TooltipString', 'Pan in X direction only', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

pan_y_image = imread('nia_pan_y_icon.png');
pan_y_ctrl = uitoggletool(tbar, ...
    'CData', pan_y_image, ...
    'TooltipString', 'Pan in y direction only', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel');

if ~isempty(live_info)
    bind_image = imread('nia_bind_icon.png');
    bind_ctrl = uitoggletool(tbar, ...
        'CData', bind_image, ...
        'TooltipString', 'Enable/disable live updates', ...
        'HandleVisibility', 'off', ...
        'Interruptible', 'off', ...
        'BusyAction', 'cancel');
else
    bind_ctrl = [];
end

% Create the plot axes
udata.axes_list = cell(1, size(data, 3));
udata.plots_list = cell(1, size(data, 3));

for idx=1:length(data)
    h = axes('Parent', fig);
    set(h, 'Units', 'Pixels');

    if idx ~= 1
        set(h, 'XTickLabel', []);
        xlabel('');
    else
        xlabel(x_label);  
    end
    
    ylabel(y_label);
    
    udata.axes_list{idx} = h;
end

% Push variables into figure user data
udata.zoom_xy_ctrl = zoom_xy_ctrl;
udata.zoom_x_ctrl = zoom_x_ctrl;
udata.zoom_y_ctrl = zoom_y_ctrl;
udata.pan_xy_ctrl = pan_xy_ctrl;
udata.pan_x_ctrl = pan_x_ctrl;
udata.pan_y_ctrl = pan_y_ctrl;
udata.bind_ctrl = bind_ctrl;
udata.zoom_h = zoom(fig);
udata.pan_h = pan(fig);
udata.live_info = live_info;
set(fig, 'UserData', udata);

% Remove the pan mode context menu
pan_context_menu = uicontextmenu;
set(udata.pan_h, 'UIContextMenu', pan_context_menu);

% Calculate the axes positions
updatePositions(fig, udata.axes_list);

% Push data to the plots
for idx=1:length(data)
    if isempty(colors)
        clr = [0.2, 0.2, 0.2];
    else
        clr = colors{idx};
    end
    
    dset = data{idx};
    
    udata.plots_list{idx} = plot(...
        udata.axes_list{idx}, dset(1,:), dset(2,:));
        
    set(udata.plots_list{idx}, ...
        'LineWidth', 2, ...
        'Color', clr);
    
    ax = udata.axes_list{idx};
    
    if idx ~= 1
        set(ax, 'XTickLabel', []);
        xlabel(ax, '');
    else
        xlabel(ax, x_label);  
    end
    
    ylabel(ax, y_label);
end

% Link the axes
linkaxes([udata.axes_list{:}], 'x');

% Store new copy of user data
set(fig, 'UserData', udata);

% Set callbacks
set(fig, ...
    'ResizeFcn', @(src, evt) figResizeCallback(fig, src, evt), ...
    'DeleteFcn', @(src, evt) figDeleteCallback(fig, src, evt));

set(zoom_orig_ctrl, ...
    'ClickedCallback', @(src, evt) zoomOrigCallback(fig, src, evt));

set(zoom_xy_ctrl, ...
    'OnCallback', @(src, evt) zoomXYOnCallback(fig, src, evt), ...
    'OffCallback', @(src, evt) zoomCtrlOffCallback(fig, src, evt));

set(zoom_x_ctrl, ...
    'OnCallback', @(src, evt) zoomXOnCallback(fig, src, evt), ...
    'OffCallback', @(src, evt) zoomCtrlOffCallback(fig, src, evt));

set(zoom_y_ctrl, ...
    'OnCallback', @(src, evt) zoomYOnCallback(fig, src, evt), ...
    'OffCallback', @(src, evt) zoomCtrlOffCallback(fig, src, evt));

set(pan_xy_ctrl, ...
    'OnCallback', @(src, evt) panXYOnCallback(fig, src, evt), ...
    'OffCallback', @(src, evt) panCtrlOffCallback(fig, src, evt));

set(pan_x_ctrl, ...
    'OnCallback', @(src, evt) panXOnCallback(fig, src, evt), ...
    'OffCallback', @(src, evt) panCtrlOffCallback(fig, src, evt));

set(pan_y_ctrl, ...
    'OnCallback', @(src, evt) panYOnCallback(fig, src, evt), ...
    'OffCallback', @(src, evt) panCtrlOffCallback(fig, src, evt));

if ~isempty(bind_ctrl)
    set(bind_ctrl, ...
        'OnCallback', @(src, evt) bindOnCallback(fig, src, evt), ...
        'OffCallback', @(src, evt) bindOffCallback(fig, src, evt));
end

end

function updatePositions(fig, axes_list)
% This function recalculates the positions of each of the axes
% based on the current dimensions of the main figure.

cfg_vert_spacing = 15;
cfg_horiz_margins = 10;
cfg_vert_margins = 10;

if isempty(axes_list)
    return;
end

fig_props = get(fig);
fig_pos = fig_props.Position;
fig_width = fig_pos(3);
fig_height = fig_pos(4);

first_inset = get(axes_list{1}, 'TightInset');

avail_width = fig_width - first_inset(1) - first_inset(3) ...
    - 2*cfg_horiz_margins;
avail_height = fig_height - first_inset(2) - 2*cfg_vert_margins ...
    - cfg_vert_spacing * (length(axes_list)-1) ...
    - first_inset(4)*length(axes_list);

plot_width = avail_width;
plot_height = avail_height / length(axes_list);

plot_width = max([1, plot_width]);
plot_height = max([1, plot_height]);

accum_ypos = first_inset(2) + cfg_vert_margins;
for idx=1:length(axes_list)
    set(axes_list{idx}, 'Position', [first_inset(1) + cfg_horiz_margins, ...
        accum_ypos, plot_width, plot_height]);
    accum_ypos = accum_ypos + cfg_vert_spacing + plot_height + ...
        first_inset(4);
end

end

function figResizeCallback(fig, ~, ~)
% This function is invoked when the figure is resized.

udata = get(fig, 'UserData');

updatePositions(fig, udata.axes_list)
end

function figDeleteCallback(fig, ~, ~)
% This function is invoked when the figure is deleted.

removeAllBindings(fig);
end

function zoomOrigCallback(fig, ~, ~)
% This function is invoked when the zooom original button is clicked

zoom(fig, 'out');
end


function zoomCtrlOffCallback(fig, ~, ~)
% This function is invoked when any zoom button is clicked off.

udata = get(fig, 'UserData');

set(udata.pan_h, 'Enable', 'off');
set(udata.zoom_h, 'Enable', 'off');

end

function zoomXYOnCallback(fig, ~, ~)
% This function is invoked when the zoom xy button is clicked on.

udata = get(fig, 'UserData');

set(udata.pan_h, 'Enable', 'off');

set(udata.zoom_h, ...
    'Enable', 'on', ...
    'Motion', 'both', ...
    'RightClickAction', 'InverseZoom');

% These do not invoke callback for respective controls, no
% idea why, matlab documentation is pretty painfully bad
set(udata.zoom_x_ctrl, 'State', 'off');
set(udata.zoom_y_ctrl, 'State', 'off');
set(udata.pan_xy_ctrl, 'State', 'off');
set(udata.pan_x_ctrl, 'State', 'off');
set(udata.pan_y_ctrl, 'State', 'off');

end

function zoomXOnCallback(fig, ~, ~)
% This function is invoked when the zoom x button is clicked on.

udata = get(fig, 'UserData');

set(udata.pan_h, 'Enable', 'off');

set(udata.zoom_h, ...
    'Enable', 'on', ...
    'Motion', 'horizontal', ...
    'RightClickAction', 'InverseZoom');

% These do not invoke callback for respective controls, no
% idea why, matlab documentation is pretty painfully bad
set(udata.zoom_xy_ctrl, 'State', 'off');
set(udata.zoom_y_ctrl, 'State', 'off');
set(udata.pan_xy_ctrl, 'State', 'off');
set(udata.pan_x_ctrl, 'State', 'off');
set(udata.pan_y_ctrl, 'State', 'off');

end

function zoomYOnCallback(fig, ~, ~)
% This function is invoked when the zoom y button is clicked on.

udata = get(fig, 'UserData');

set(udata.pan_h, 'Enable', 'off');

set(udata.zoom_h, ...
    'Enable', 'on', ...
    'Motion', 'vertical', ...
    'RightClickAction', 'InverseZoom');

% These do not invoke callback for respective controls, no
% idea why, matlab documentation is pretty painfully bad
set(udata.zoom_xy_ctrl, 'State', 'off');
set(udata.zoom_x_ctrl, 'State', 'off');
set(udata.pan_xy_ctrl, 'State', 'off');
set(udata.pan_x_ctrl, 'State', 'off');
set(udata.pan_y_ctrl, 'State', 'off');

end

function panCtrlOffCallback(fig, ~, ~)
% This function is invoked when any pan button is clicked off.

udata = get(fig, 'UserData');

set(udata.pan_h, 'Enable', 'off');
set(udata.zoom_h, 'Enable', 'off');

end

function panXYOnCallback(fig, ~, ~)
% This function is invoked when the pan xy button is clicked on.

udata = get(fig, 'UserData');

set(udata.zoom_h, 'Enable', 'off');

set(udata.pan_h, ...
    'Enable', 'on', ...
    'Motion', 'both');

% These do not invoke callback for respective controls, no
% idea why, matlab documentation is pretty painfully bad
set(udata.zoom_xy_ctrl, 'State', 'off');
set(udata.zoom_x_ctrl, 'State', 'off');
set(udata.zoom_y_ctrl, 'State', 'off');
set(udata.pan_x_ctrl, 'State', 'off');
set(udata.pan_y_ctrl, 'State', 'off');

end

function panXOnCallback(fig, ~, ~)
% This function is invoked when the pan x button is clicked on.

udata = get(fig, 'UserData');

set(udata.zoom_h, 'Enable', 'off');

set(udata.pan_h, ...
    'Enable', 'on', ...
    'Motion', 'horizontal');

% These do not invoke callback for respective controls, no
% idea why, matlab documentation is pretty painfully bad
set(udata.zoom_xy_ctrl, 'State', 'off');
set(udata.zoom_x_ctrl, 'State', 'off');
set(udata.zoom_y_ctrl, 'State', 'off');
set(udata.pan_xy_ctrl, 'State', 'off');
set(udata.pan_y_ctrl, 'State', 'off');

end

function panYOnCallback(fig, ~, ~)
% This function is invoked when the pan y button is clicked on.

udata = get(fig, 'UserData');

set(udata.zoom_h, 'Enable', 'off');

set(udata.pan_h, ...
    'Enable', 'on', ...
    'Motion', 'vertical');

% These do not invoke callback for respective controls, no
% idea why, matlab documentation is pretty painfully bad
set(udata.zoom_xy_ctrl, 'State', 'off');
set(udata.zoom_x_ctrl, 'State', 'off');
set(udata.zoom_y_ctrl, 'State', 'off');
set(udata.pan_xy_ctrl, 'State', 'off');
set(udata.pan_x_ctrl, 'State', 'off');

end

function bindOnCallback(fig, ~, ~)
% This function is invoked when binding (live update) are enabled.

udata = get(fig, 'UserData');

if isempty(udata.live_info)
    return;
end

if ~ishghandle(udata.live_info.creator)
    return;
end

creator_udata = get(udata.live_info.creator, 'UserData');

for axes_idx=1:length(udata.axes_list)
    creator_udata.addListenerToROI(udata.live_info.uuids{axes_idx}, ...
        udata.plots_list{axes_idx});
end

end

function removeAllBindings(fig)
% This function removes all bindings (live updates).

udata = get(fig, 'UserData');

if isempty(udata.live_info)
    return;
end

if ~ishghandle(udata.live_info.creator)
    return;
end

creator_udata = get(udata.live_info.creator, 'UserData');

for axes_idx=1:length(udata.axes_list)
    creator_udata.removeListenerFromROI(udata.live_info.uuids{axes_idx}, ...
        udata.plots_list{axes_idx});
end

end

function bindOffCallback(fig, ~, ~)
% This function is invoked when binding (live update) are disabled.

removeAllBindings(fig);
end



