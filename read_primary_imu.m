function varargout = read_primary_imu(filepath, max_records, mode)
% READ_PRIMARY_IMU Reads binary IMU data from file line-by-line.
%
% Interprets binary data according to epost.pdf:
%   struct Generic_IMU_record {
%       double Tow;    // Time of Week          [ s ]
%       int    dw[3];  // Angular increments    [pulse count]
%       int    dv[3];  // Velocity increments   [pulse count]
%   };                 // Size: 32 B
%
% Conversion factors:
%   - Gyro:          206264806.247 pulses/rad  -> dw [rad]
%   - Accelerometer: 20000000.0    pulses/(m/s) -> dv [m/s]
%
% Usage:
%   data = read_primary_imu()
%   data = read_primary_imu(filepath)
%   data = read_primary_imu(filepath, max_records)
%   data = read_primary_imu(filepath, max_records, mode)
%   [tow, dw, dv] = read_primary_imu(...)
%   [tow, dw, dv, raw] = read_primary_imu(...)
%   [data, raw] = read_primary_imu(...)
%
% Inputs:
%   filepath    - Path to binary file (default: 'data/primary_imu.dat')
%   max_records - Maximum number of records to read (default: inf, all)
%   mode        - Reading mode:
%                   'line' (default) - reads record by record (line by line)
%                   'fast' / 'vectorized' - reads in bulk for maximum speed
%
% Outputs:
%   When 1 output requested:
%     data - Struct with fields:
%              .tow   - [N x 1 double] Time of Week [s]
%              .dw    - [N x 3 double] Angular increments [rad]
%              .dv    - [N x 3 double] Velocity increments [m/s]
%              .raw   - Struct with raw integer pulse counts (.dw, .dv)
%              .dt    - Nominal / average sample interval [s]
%              .fs    - Estimated sample frequency [Hz]
%              .units - Description of physical units
%
%   When 2 outputs requested:
%     [data, raw]
%
%   When 3 or more outputs requested:
%     [tow, dw, dv, raw] or [tow, dw, dv, raw_dw, raw_dv]

    % Default arguments
    if nargin < 1 || isempty(filepath)
        filepath = 'data/primary_imu.dat';
    end

    % Resolve file path if relative path is not found in current folder
    if ~exist(filepath, 'file')
        script_dir = fileparts(mfilename('fullpath'));
        alt_path = fullfile(script_dir, filepath);
        if exist(alt_path, 'file')
            filepath = alt_path;
        elseif exist(fullfile(script_dir, 'data', 'primary_imu.dat'), 'file')
            filepath = fullfile(script_dir, 'data', 'primary_imu.dat');
        else
            error('read_primary_imu:FileNotFound', 'File not found: %s', filepath);
        end
    end

    if nargin < 2 || isempty(max_records)
        max_records = inf;
    end

    if nargin < 3 || isempty(mode)
        mode = 'line';
    end

    % Scale factors from epost.pdf
    SCALE_GYRO = 206264806.247; % pulses / rad
    SCALE_ACC  = 20000000.0;     % pulses / (m/s)
    RECORD_SIZE = 32;            % 8 (double) + 3*4 (int32) + 3*4 (int32) bytes

    % Open file
    fid = fopen(filepath, 'rb');
    if fid == -1
        error('read_primary_imu:FileOpenFailed', 'Could not open file: %s', filepath);
    end
    cleaner = onCleanup(@() fclose(fid));

    % Determine total number of records available in file
    fseek(fid, 0, 'eof');
    file_bytes = ftell(fid);
    fseek(fid, 0, 'bof');

    total_records = floor(file_bytes / RECORD_SIZE);
    n_to_read = min(total_records, max_records);

    if strcmpi(mode, 'fast') || strcmpi(mode, 'vectorized')
        % Bulk read mode
        raw_bytes = fread(fid, n_to_read * RECORD_SIZE, '*uint8');
        actual_records = floor(length(raw_bytes) / RECORD_SIZE);

        raw_mat = reshape(raw_bytes(1:(actual_records * RECORD_SIZE)), RECORD_SIZE, actual_records);
        tow = typecast(reshape(raw_mat(1:8, :), 1, []), 'double')';
        raw_dw = reshape(typecast(reshape(raw_mat(9:20, :), 1, []), 'int32'), 3, actual_records)';
        raw_dv = reshape(typecast(reshape(raw_mat(21:32, :), 1, []), 'int32'), 3, actual_records)';
    else
        % Line-by-line reading mode (reading record by record)
        tow = zeros(n_to_read, 1);
        raw_dw = zeros(n_to_read, 3, 'int32');
        raw_dv = zeros(n_to_read, 3, 'int32');

        actual_records = 0;
        for k = 1:n_to_read
            t = fread(fid, 1, 'double');
            if isempty(t)
                break;
            end
            w = fread(fid, 3, 'int32');
            v = fread(fid, 3, 'int32');

            if length(w) < 3 || length(v) < 3
                break;
            end

            actual_records = actual_records + 1;
            tow(actual_records) = t;
            raw_dw(actual_records, :) = w;
            raw_dv(actual_records, :) = v;
        end

        % Truncate if fewer records than expected
        if actual_records < n_to_read
            tow = tow(1:actual_records);
            raw_dw = raw_dw(1:actual_records, :);
            raw_dv = raw_dv(1:actual_records, :);
        end
    end

    % Convert raw pulses to physical quantities
    dw = double(raw_dw) / SCALE_GYRO;
    dv = double(raw_dv) / SCALE_ACC;

    % Calculate sampling information
    if actual_records > 1
        dt = (tow(end) - tow(1)) / (actual_records - 1);
        fs = 1 / dt;
    else
        dt = NaN;
        fs = NaN;
    end

    raw_struct = struct('dw', raw_dw, 'dv', raw_dv);

    data = struct();
    data.tow = tow;
    data.dw = dw;
    data.dv = dv;
    data.raw = raw_struct;
    data.dt = dt;
    data.fs = fs;
    data.units = struct(...
        'tow', 'seconds', ...
        'dw',  'radians', ...
        'dv',  'm/s', ...
        'raw_dw', 'pulses', ...
        'raw_dv', 'pulses' ...
    );

    % Assign requested outputs
    if nargout <= 1
        varargout{1} = data;
    elseif nargout == 2
        varargout{1} = data;
        varargout{2} = raw_struct;
    elseif nargout == 3
        varargout{1} = tow;
        varargout{2} = dw;
        varargout{3} = dv;
    elseif nargout == 4
        varargout{1} = tow;
        varargout{2} = dw;
        varargout{3} = dv;
        varargout{4} = raw_struct;
    elseif nargout >= 5
        varargout{1} = tow;
        varargout{2} = dw;
        varargout{3} = dv;
        varargout{4} = raw_dw;
        varargout{5} = raw_dv;
    end
end
