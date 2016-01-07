classdef nia_movie < handle
    %NIA_MOVIE Movie object
    %   The nia_movie class stores information on a movie. Movies are
    %   broken into a series of "slices." Each slice has a timestamp,
    %   a position vector, and a list of associated channel data entries.
    %   Each channel data entry has a channel identifier and a chunk of
    %   image data. Note that nia_movie is a handle-based object, and
    %   so it follows handle semantics and not value semantics.
    
    properties (Access=protected, Hidden=true)
        
        % This property is a structure array of slices. Each slice has the
        % fields 'time', 'channels', and 'pos'. The time field is a scalar
        % double. The 'channels' field is a structure array with the fields
        % 'ch' and 'image', where the 'ch' field is a channel identifier,
        % and the 'image' field is a image blob. The 'pos' field of slices
        % is the position vector. The position and channel identifier must
        % have positive integer values.
        slices;
        
        % This property is [1xN] array of channel identifiers. Every
        % channel identifier used in the movie is kept here in a sorted
        % list. This property is optional and may be empty.
        ch_list;
        
        % This property is a [1xM] vector where is M is the number of
        % position vectors that contains the length of each position
        % vector. All slices must have the same length for each position
        % vector.
        pos_lens;
        
        % This property a lookup table used to locate first slice in the
        % movie with a given position. The format of this property is
        % described in loadPosCache(). This property is optional and may be
        % empty.
        pos_lu;
        
        % This property is a [1xM] cell array of [2xN] arrays, where M is
        % the number of position vectors and N is the length of each
        % position vector. Each column specifies the inclusive range for
        % values in a given position dimension, with the first row storing
        % the minimum value and the second row storing the maximum value.
        % This property is optional and may be empty.
        pos_ranges;
        
        % This property structure array that specifies the histogram
        % information for each channel of the movie. It has a number of
        % elements equal to the number of channels (the length of ch_list)
        % each with the fields 'min', 'max', and 'dist'. The min field is
        % scalar storing the minimum intensity value, the max field is a
        % scalar storing the maximum intensity value, and dist is a [2xN]
        % array specifying the histogram count for a series of evenly
        % spaced bins. Each column of dist describes a bin, with the center
        % position of the bin in the first row, and the count in the second
        % row. This property is optional and may be empty.
        hist_info;
        
        % This property stores the interface version for the movie. The
        % importStruct() function will only accept structure movies with a
        % lower version number.  This number should be incremented every
        % time there is a backwards incompatible change.
        api_ver = 1;
        
        % This property stores a boolean value indicating if the movie is
        % operating in "safe mode," where all input are arguments are
        % carefully verified.
        safe_mode = true;
    end
    
    methods (Access=public)
        function obj = nia_movie()
            % Create an empty movie object
            
            obj.slices = [];
            obj.ch_list = [];
            obj.pos_lens = [];
            obj.pos_lu = {};
            obj.pos_ranges = {};
            obj.hist_info = [];
            % leave api_ver at default above
            % leave safe_mode at default above
        end
        
        out = copy(obj);
    end
    
    methods (Access=public)
        loadFlatMovie(obj, fmovie, freq)
        loadPrairieTSeriesXML(obj, folder, fname)
        
        loadChCache(obj)
        loadPosCache(obj)
        loadHistCache(obj)
        
        invalidatePosCache(obj);
        
        state = isEmpty(obj)
        num = getNumSlices(obj)
        pos_len = getPositionLengths(obj)
        ch_list = getChannelList(obj)
        
        slice_info = getSliceInfo(obj, slice_idx)
        setSliceInfo(obj, slice_idx, slice_info, force_update)
        
        slice_idx = appendSlice(obj, slice_info)
        removeSlice(obj, slice_idx)
        
        ch_data = getChannelData(obj, slice_idx, ch_idx)
        ch_data = getChannelDataById(obj, slice_idx, ch_id)
        
        setChannelData(obj, slice_idx, ch_idx, ch_data)
        setChannelDataById(obj, slice_idx, ch_id, ch_data)
        
        ch_idx = addChannelData(obj, slice_idx, ch_data)
        
        removeChannelData(obj, slice_idx, ch_idx)
        removeChannelDataById(obj, slice_idx, ch_id)
        
        pos_ranges = getPositionRanges(obj)
        slice_idx = getSliceIndex(obj, pos_vec, pos)
        hist_info = getHistInfo(obj)
        
        [out, msg] = isValidRange(obj, rnge)
        
        process(obj, rnge, slice_transform, ch_transform)
        append(obj, mov2, rnge, transform)
        [udata, stop] = analyze(obj, rnge, transform, udata)
        
        [accum, count] = accumulate(obj, rnge, transform)
        [flat, info] = flatten(obj, rnge, transform)
        
        out = scan(obj, pos_vec, home_pos, scan_extents, channel_id, ...
            transform, udata);
        out = scanToLast(obj, pos_vec, home_pos, scan_extents);
        out = scanROI(obj, pos_vec, home_pos, scan_pos, channel_id, mask)
        
        sdata = exportStruct(obj)
        importStruct(obj, sdata)
    end
    
    methods (Access=protected, Hidden=true)
        rnge = normalizeRange(obj, rnge)
        allowed = inSliceRange(obj, slice_idx, rnge)
        allowed = inChannelRange(obj, slice_idx, ch_idx, rnge)
    end
end
