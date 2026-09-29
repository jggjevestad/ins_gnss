function [rec, is_eof] = read_imu_record(fid)
% READ_IMU_RECORD Reads a single line/record from an open binary IMU file.
%
% Interpreted according to epost.pdf:
%   struct Generic_IMU_record {
%       double Tow;    // Time of Week          [ s ]
%       int    dw[3];  // Angular increments    [pulse count]
%       int    dv[3];  // Velocity increments   [pulse count]
%   };                 // Total size: 32 bytes
%
% Conversion factors:
%   - Gyro:          206264806.247 pulses/rad
%   - Accelerometer: 20000000.0 pulses/(m/s)
%
% Syntax:
%   [rec, is_eof] = read_imu_record(fid)
%
% Inputs:
%   fid    - File identifier from fopen(..., 'rb')
%
% Outputs:
%   rec    - Struct with fields:
%              .tow    - Time of Week [s] (double)
%              .dw     - Angular increments [rad] (1x3 double)
%              .dv     - Velocity increments [m/s] (1x3 double)
%              .raw_dw - Raw angular pulses (1x3 int32)
%              .raw_dv - Raw velocity pulses (1x3 int32)
%   is_eof - Boolean flag indicating if end of file was reached

    SCALE_GYRO = 206264806.247; % pulses / rad
    SCALE_ACC  = 20000000.0;     % pulses / (m/s)

    rec = struct('tow', [], 'dw', [], 'dv', [], 'raw_dw', [], 'raw_dv', []);
    is_eof = false;

    tow = fread(fid, 1, 'double');
    if isempty(tow)
        is_eof = true;
        return;
    end

    raw_dw = fread(fid, 3, 'int32');
    raw_dv = fread(fid, 3, 'int32');

    if length(raw_dw) < 3 || length(raw_dv) < 3
        is_eof = true;
        return;
    end

    rec.tow = tow;
    rec.raw_dw = raw_dw(:)';
    rec.raw_dv = raw_dv(:)';
    rec.dw = double(raw_dw(:)') / SCALE_GYRO;
    rec.dv = double(raw_dv(:)') / SCALE_ACC;
end
