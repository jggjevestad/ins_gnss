% Main script for INS/GNSS processing
% Demonstrates reading and inspecting primary IMU binary data

clear; clc;

imu_file = fullfile('data', 'primary_imu.dat');

% Read first 10,000 records line-by-line (or omit 10000 to read all records)
fprintf('Reading IMU data from: %s\n', imu_file);
data = read_primary_imu(imu_file, 10000);

fprintf('\n--- IMU Data Summary ---\n');
fprintf('Records read:      %d\n', length(data.tow));
fprintf('Start GPS ToW:     %.3f s\n', data.tow(1));
fprintf('End GPS ToW:       %.3f s\n', data.tow(end));
fprintf('Duration:          %.3f s\n', data.tow(end) - data.tow(1));
fprintf('Sample rate:       %.2f Hz (dt = %.5f s)\n', data.fs, data.dt);
fprintf('Mean acceleration: [%.3f, %.3f, %.3f] m/s^2\n', mean(data.dv) / data.dt);
fprintf('Mean angular rate: [%.4f, %.4f, %.4f] deg/s\n', rad2deg(mean(data.dw) / data.dt));