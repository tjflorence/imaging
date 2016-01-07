function loadPrairieTSeriesXML(obj, folder, fname)
%LOADPRAIRIETSERIESXML Read movie from priaire xml file
%   loadPrairieTSeriesXML(folder, fname) loads the XML file in the folder
%   'folder' with the filename 'fname.' Any existing movie contents are
%   destroyed. The XML file should describe a series of sequences, each
%   containing a series of channels. If the XML describes a TSeries with
%   multiple sequences, then the first element of the position vector for
%   each slice will contain the index of the sequence and the second
%   element of the position vector will contain the index of the frame in
%   the sequence. If the XML descrbies a TSeries with only a single
%   sequence, then the position vector contains only a single element
%   describing the index of the frame in sequence.
%   
%   Performance note: MATLAB 2014a has a bug that causes it to consume
%   large amounts of CPU time if the current directory contains many files
%   (possibly related to the misguided inclusion of the current directory
%   in the path). As such, it is best to use the 'folder' argument to
%   specify the folder containing the TIFF stack rather than changing
%   directories into the folder containing the TIFF stack.

% Check input parameters
if ~nia_isString(folder)
    error 'The argument ''folder'' must be a string';
end

if ~nia_isString(fname)
    error 'The argument ''fname'' must be string'
end

% Destroy the movie
obj.slices = [];
obj.ch_list = [];
obj.pos_lens = [];
obj.pos_lu = {};
obj.pos_ranges = {};
obj.hist_info = [];

% Find basename of file
if ~strcmp(fname(end-3:end), '.xml')
    error 'Argument fname has unrecognized suffix'
end

fname_base = fname(1:end-4);

% Find name of cleaned version
fname_full = [folder, filesep, fname];
fname_clean = [folder, filesep, fname_base, '_trimmed.xml'];

% If unclean version doesn't exist, then error
if ~exist(fname_full, 'file')
    error('The file ''%s'' does not exist', fname_full);
end

% If clean version does not exist, create it
if ~exist(fname_clean, 'file')
    xslt(fname_full, 'nia_trim_prairie_tseries_xml.xml', fname_clean);
end

% Read in the XML file
try
    root = xmlread(fname_clean);
catch
    error 'Unable to parse XML file';
end

% Search for the PVScan element
if ~root.hasChildNodes
    return;
end

top_children = root.getChildNodes;
for idx=1:top_children.getLength
    if strcmp(top_children.item(idx-1).getNodeName, 'PVScan')
        pvscan_root = top_children.item(idx-1);
        break;
    end
end

% Get nodes of the sequence element
if ~pvscan_root.hasChildNodes
    return;
end

pvscan_children = pvscan_root.getChildNodes;

% Count the total number of slices
max_total_slices = 0;

for pvscan_idx=1:pvscan_children.getLength
    % Skip non-Sequence elements
    if ~strcmp(pvscan_children.item(pvscan_idx-1).getNodeName, 'Sequence')
        continue;
    end
    
    seq_root = pvscan_children.item(pvscan_idx-1);
    
    % Get nodes of the sequence element
    if ~seq_root.hasChildNodes
        continue;
    end
    
    seq_children = seq_root.getChildNodes;
    max_total_slices = max_total_slices + seq_children.getLength;
end
    
% Initialize channel list and indices
ch_list = [];
seq_cur_idx = 1;
frame_cur_idx = 1;
slice_cur_idx = 1;

% Pre-allocate the slices array (it is too large
% but we will trim at the end)
slices(max_total_slices).time = [];
slices(max_total_slices).channels = [];
slices(max_total_slices).pos = {};

% Traverse list of children
for pvscan_idx=1:pvscan_children.getLength
    
    % Skip non-Sequence elements
    if ~strcmp(pvscan_children.item(pvscan_idx-1).getNodeName, 'Sequence')
        continue;
    end
    
    seq_root = pvscan_children.item(pvscan_idx-1);
    
    % Get nodes of the sequence element
    if ~seq_root.hasChildNodes
        seq_cur_idx = seq_cur_idx + 1;
        continue;
    end
    
    seq_children = seq_root.getChildNodes;
    
    % Traverse the list of children
    for seq_idx=1:seq_children.getLength
        
        % Skip non-Frame elements
        if ~strcmp(seq_children.item(seq_idx-1).getNodeName, 'Frame')
            continue;
        end
        
        frame_root = seq_children.item(seq_idx-1);
        
        % Parse frame attributes
        attrib_list = frame_root.getAttributes;
        
        t_val = [];
        for attrib_idx=1:attrib_list.getLength
            attrib = attrib_list.item(attrib_idx-1);
            
            if strcmp(attrib.getName, 'relativeTime')
                t_str = char(attrib.getValue);
                t_val = str2double(t_str);
                if isempty(t_val)
                    error 'Unable to parse relativeTime value';
                end
            end
        end
        
        if isempty(t_val)
            error 'The relativeTime attribute is missing';
        end
        
        slices(slice_cur_idx).time = t_val;
        slices(slice_cur_idx).pos{1} = [seq_cur_idx, frame_cur_idx];
        
        % Get nodes of the frame element
        if ~frame_root.hasChildNodes
            error 'Frames element is missing File element';
        end
        
        frame_children = frame_root.getChildNodes;
        
        % Preallocate channels array (it is too large
        % but we will trim it at the end)
        clear channels;
        channels(frame_children.getLength).ch = []; %#ok<AGROW>
        channels(frame_children.getLength).image = []; %#ok<AGROW>
        channel_cur_idx = 1;
        
        for frame_idx=1:frame_children.getLength
            
            % Skip non-File elements
            if ~strcmp(frame_children.item(frame_idx-1).getNodeName, 'File')
                continue;
            end
            
            file_root = frame_children.item(frame_idx-1);
            attrib_list = file_root.getAttributes;
            
            ch_val = [];
            fname_val = [];
            for attrib_idx=1:attrib_list.getLength
                attrib = attrib_list.item(attrib_idx-1);
                
                if strcmp(attrib.getName, 'channel')
                    if ~isempty(ch_val)
                        error 'File element contains redundant channel attribute';
                    end
                    
                    a_str = char(attrib.getValue);
                    ch_val = str2double(a_str);
                    
                    if isempty(ch_val)
                        error 'Unable to parse channel attribute';
                    end
                elseif strcmp(attrib.getName, 'filename')
                    if ~isempty(fname_val)
                        error 'File element contains redundant filename attribute';
                    end
                    
                    fname_val = char(attrib.getValue);
                end
            end
            
            if isempty(ch_val)
                error 'File element is missing channel attribute';
            end
            
            if isempty(fname_val)
                error 'File element is missing filename attribute';
            end
            
            fname_val_full = [folder, filesep, fname_val];
            
            channels(channel_cur_idx).ch = ch_val; %#ok<AGROW>
            channels(channel_cur_idx).image = double(imread(fname_val_full)); %#ok<AGROW>
            channel_cur_idx = channel_cur_idx + 1;
        end
        
        % Trim the channels array
        if channel_cur_idx <= 1
            channels = [];
        else
            clear tmp;
            tmp(channel_cur_idx-1).ch = []; %#ok<AGROW>
            tmp(channel_cur_idx-1).image = []; %#ok<AGROW>
            tmp(1:channel_cur_idx-1) = channels(1:channel_cur_idx-1);
            channels = tmp;
        end
        
        % Check that channels list is unique
        if ~isempty(channels)
            frame_ch_list = [channels(:).ch];
            if length(frame_ch_list) < length(unique(frame_ch_list))
                error 'Frame element contains redundant channel entries';
            end
            
            ch_list = unique([ch_list, frame_ch_list]);
        end
        
        slices(slice_cur_idx).channels = channels;
        slice_cur_idx = slice_cur_idx + 1;
        frame_cur_idx = frame_cur_idx + 1;
    end
    
    seq_cur_idx = seq_cur_idx + 1;
    frame_cur_idx = 1;
end

% Trim the slices array
if slice_cur_idx <= 1
    slices = [];
else
    slices(slice_cur_idx:end) = [];
end

% If only one sequence, then remove position indices
if seq_cur_idx <= 2
    for slice_idx=1:length(slices)
        slices(slice_idx).pos{1} = slices(slice_idx).pos{1}(2); %#ok<AGROW>
    end
    
    pos_lens_entry = 1;
else
    pos_lens_entry = 2;
end

obj.slices = slices;
obj.ch_list = ch_list;
obj.pos_lens = pos_lens_entry;

end
