% INS GNSS example

% Notation
% ECEF: Earth-Centered, Earth-Fixed frame
% NED: North-East-Down frame
% Platform: platform/IMU frame
% Body: body frame
% Camera: camera frame
% Map: map frame
% 
% Rotation matrices
% Ce_g: ECEF frame to local level geodetic (NED) frame
% Cg_p: local level geodetic (NED) frame to platform/IMU frame (roll/pitch/yaw)
% Cs_p: sensor to platform/IMU frame
% Cp_b: platform/IMU frame to body frame (for now Cp_b = identity)
% Cp_c: platform/IMU frame to camera frame (not used in this simple example)
% Cm_c: map to camera frame (not used in this simple example)

% Main function
function main
    % System parameters
    filename = 'data/primary_imu.dat';
    N = 1000;

    % Initial parameters from first line of data/traj.txt (epoch 198300.000259 s)
    t0   = 198300.000259;                              % Initial GPS Time of Week [s]
    x0_e = [3208319.5845; 633237.6578; 5458116.3716];  % Initial position in ECEF [m]
    v0_e = [-54.8084; -4.6903; 32.4978];               % Initial velocity in ECEF [m/s]
    Cg_p = [ 0.994637,  0.099168, -0.029369;           % Local geodetic (NED) to platform frame
            -0.102838,  0.918074, -0.382837;
            -0.011003,  0.383804,  0.923349];

    % Read IMU observations from binary file starting at epoch t0
    [dw, dv, t] = read_imu(filename, N, t0);

    fprintf('Initial epoch ToW: %.6f s\n', t(1));
    fprintf('Read %d IMU records.\n', size(dw, 1));
    fprintf('dw (rad): [%.3e, %.3e, %.3e]\n', dw(1, :));
    fprintf('dv (m/s): [%.3e, %.3e, %.3e]\n', dv(1, :));
end

% Read IMU observations from binary file starting from epoch t0
function [dw, dv, t] = read_imu(filename, N, t0)
    fid = fopen(filename, 'rb');

    % Seek to epoch t0 if specified
    if nargin >= 3 && ~isempty(t0)
        t1 = fread(fid, 1, 'double');
        idx = max(0, round((t0 - t1) / 0.005));
        fseek(fid, idx * 32, 'bof');
        t_curr = fread(fid, 1, 'double');
        while abs(t_curr - t0) > 1e-4
            step = round((t0 - t_curr) / 0.005);
            if step == 0, break; end
            idx = max(0, idx + step);
            fseek(fid, idx * 32, 'bof');
            t_curr = fread(fid, 1, 'double');
        end
        fseek(fid, idx * 32, 'bof');
    end

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

% Skew-symmetric matrix from vector
function crossm = skew(x)
    crossm = [0, -x(3), x(2);
              x(3), 0, -x(1);
              -x(2), x(1), 0];
end