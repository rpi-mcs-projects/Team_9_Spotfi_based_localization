function data_struct = load_h5_structure(filepath, num_data, print_struct_flag)
   
    if nargin < 2
        error('load_h5_structure:MissingNumData', ...
            ['second argument "num_data" is REQUIRED to limit the data length to load']);
    end

    if ~isscalar(num_data) || ~isnumeric(num_data) || ~isfinite(num_data) ...
            || num_data <= 0 || floor(num_data) ~= num_data
        error('load_h5_structure:InvalidNumData', ...
            '"num_data" must be a finite positive integer.');
    end
    
    if print_struct_flag
        info = h5info(filepath);
        print_h5_structure(info,'');
    end

    start = [1, 1, 1, 1];
    count = [num_data, inf, inf, inf];

    real_channel = h5read(filepath, '/channels/real', start, count); % [N, F, ant, AP]
    imag_channel = h5read(filepath, '/channels/imag', start, count); % [N, F, ant, AP]

    % max N = 17276
    % max F = 234
    % max ant = 4
    % max AP = 6

    data_struct.cmplx_channel = complex(real_channel, imag_channel); % [num_data, F, ant, AP]
    data_struct.ap_locs = h5read(filepath, '/AP_locs'); % [AP, ant, 2], location of 6 APs and each of its 4 antennas
    data_struct.labels = h5read(filepath, '/labels', [1, 1], [num_data, inf]); %[num_data, 2] GT locations
    data_struct.rssi = h5read(filepath, '/rssi', [1,1], [num_data, inf]); % [num_data, AP], rssi at each AP for each data
    data_struct.ant_sep = h5read(filepath, '/opt/ANT_SEP'); % [1,1] scalar, antenna separation on APs
    data_struct.bandwidth = h5read(filepath, '/opt/BW'); % [1,1] scalar, bandwidth
    data_struct.center_freq = h5read(filepath, '/opt/CENTER_FREQ'); % [1,1] scalar, center freq
    data_struct.wifi_channel = h5read(filepath, '/opt/CHAN'); % [1,1] scalar wifi 11ac channel number used 
    data_struct.sub_carrier_idx = h5read(filepath, '/opt/SUB_IND'); % [1, F] subcarrier indices used for Wifi
    data_struct.freq_list = h5read(filepath, '/opt/FREQ'); % [1, F] subcarrier frequencies used for WiFi

end


function print_h5_structure(info, prefix)
% Recursive helper to print HDF5 structure
    for i = 1:length(info.Groups)
        fprintf('%s[GROUP] %s:\n', prefix, info.Groups(i).Name);
        print_h5_structure(info.Groups(i), [prefix '  ']);
    end
    for i = 1:length(info.Datasets)
        ds = info.Datasets(i);
        shape_str = sprintf('%dx', ds.Dataspace.Size);
        shape_str = shape_str(1:end-1);
        fprintf('%s  %s shape: [%s]  dtype: %s specific type: %s\n', prefix, ds.Name, shape_str, ds.Datatype.Class, ds.Datatype.Type);
    end
end