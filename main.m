% INS GNSS example

% Main function
function main
    % System parameters
    filename = 'data/primary_imu.dat';
    N = 1000;

    % Read IMU observations from binary file
    [dw, dv, t] = read_imu(filename, N);

    fprintf('Read %d IMU records.\n', size(dw, 1));
    fprintf('dw (rad): [%.3e, %.3e, %.3e]\n', dw(1, :));
    fprintf('dv (m/s): [%.3e, %.3e, %.3e]\n', dv(1, :));
end

% Read IMU observations from binary file
function [dw, dv, t] = read_imu(filename, N)
    fid = fopen(filename, 'rb');
    t  = zeros(N, 1);
    dw = zeros(N, 3);
    dv = zeros(N, 3);

    for k = 1:N
        t(k)     = fread(fid, 1, 'double');
        dw(k, :) = fread(fid, 3, 'int32');
        dv(k, :) = fread(fid, 3, 'int32');
    end
    fclose(fid);

    dw = dw / 206264806.247; % rad
    dv = dv / 20000000.0;     % m/s
end
