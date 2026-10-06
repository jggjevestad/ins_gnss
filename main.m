% INS GNSS example
%
% Notation
% ECEF: Earth-Centered, Earth-Fixed frame (e-frame)
% NED:  Local level geodetic frame (North-East-Down, g-frame)
% Platform/Sensor: IMU body/sensor frame (s-frame = p-frame)
%
% Direction Cosine Matrices (DCM):
% Ce_g: ECEF to local level (NED) frame
% Cg_p: NED to platform frame (roll/pitch/yaw)
% Cs_p: Sensor to platform frame
% Cp_b: Platform to body frame
% Cp_c: Platform to camera frame (sensor orientation)
% Cm_c: Map to camera frame (omega/phi/kappa)

function main
    % System parameters
    filename = 'data/primary_imu.dat';
    N = 1000;

    % Initial parameters from first line of data/traj.txt (epoch 198300.000259 s)
    t0   = 198300.000259;                              % Initial GPS Time of Week [s]
    x0_e = [3208319.5845; 633237.6578; 5458116.3716];  % Initial position in ECEF [m]
    v0_e = [-54.8084; -4.6903; 32.4978];               % Initial velocity in ECEF [m/s]
    Cg_p = [ 0.994637,  0.099168, -0.029369;           % NED to platform frame
            -0.102838,  0.918074, -0.382837;
            -0.011003,  0.383804,  0.923349];

    % Read IMU observations
    [dw, dv, t] = read_imu(filename, N, t0);
    fprintf('Initial epoch ToW: %.6f s | Records: %d\n', t(1), size(dw, 1));
    fprintf('dw (rad): [%.3e, %.3e, %.3e] | dv (m/s): [%.3e, %.3e, %.3e]\n', dw(1, :), dv(1, :));

    % Earth rotation rate vector and skew matrix in ECEF [rad/s]
    wie_e = [0; 0; 7292115e-11];
    Wie_e = skew(wie_e);

    % Initial attitude DCM: Cs_e = (Cg_p * Ce_g)'
    Cs_e = (Cg_p * ecef2ned(x0_e))';

    % State trajectories
    x_est = zeros(N, 3);
    v_est = zeros(N, 3);
    x_est(1, :) = x0_e';
    v_est(1, :) = v0_e';

    x = x0_e;
    v = v0_e;

    fprintf('\nPropagating navigation states (DCM mechanization, C.2 Eq. 59-62b)...\n');
    tic;
    for k = 1:N-1
        dt = t(k+1) - t(k);
        if dt <= 0, dt = 0.005; end

        dw_k = dw(k, :)';
        dv_k = dv(k, :)';

        % DCM attitude integration (Eq. 31, 61)
        [Cs_e_next, rho] = integrate_dcm(Cs_e, dw_k, dt, wie_e);

        % Numerical integration of position and velocity (C.2, Eq. 59-62b)
        [x, v, ~] = integrate_pos_vel(x, v, dv_k, rho, dt, Cs_e_next, Wie_e);
        Cs_e = Cs_e_next;

        x_est(k+1, :) = x';
        v_est(k+1, :) = v';
    end
    t_elapsed = toc;

    fprintf('Propagation completed: %d epochs in %.4f s (%.1f kHz)\n', ...
        N, t_elapsed, (N - 1) / t_elapsed / 1000);
    fprintf('Estimated position (end) [m]:   [%.4f, %.4f, %.4f]\n', x_est(end, :));
    fprintf('Estimated velocity (end) [m/s]: [%.4f, %.4f, %.4f]\n', v_est(end, :));

    % Benchmark against reference trajectory if available
    compare_reference(t, x_est, v_est, 'data/traj.txt');
end

% C.2 Numerical integration of position and velocity (Eq. 59-62b)
function [x, v, dv_e] = integrate_pos_vel(x, v, dv, rho, dt, Cs_e_next, Wie_e)
    % Eq. (60): Midpoint average DCM over update interval
    Cbar = Cs_e_next * (eye(3) - 0.5 * skew(rho));

    % Normal gravity vector in ECEF frame: ge = C_g^e * gamma^g
    ge = normal_gravity_ecef(x);

    % Eq. (59): Velocity increment transformed to ECEF with Coriolis and gravity
    dv_e = Cbar * dv - (2 * Wie_e * v - ge) * dt;

    % Eq. (62a-b): Position and velocity updates
    x = x + v * dt + dv_e * (dt / 2);
    v = v + dv_e;
end

% DCM attitude integration (Eq. 31, 61) with Rodrigues rotation and re-orthogonalization
function [Cs_e_next, rho] = integrate_dcm(Cs_e, dw, dt, wie_e)
    % Eq. (61): Rotation increment relative to ECEF frame
    rho = dw - Cs_e' * wie_e * dt;
    mag = norm(rho);
    sk  = skew(rho);

    % Relative DCM increment (matrix exponential via Rodrigues formula)
    if mag > 1e-12
        dC = eye(3) + (sin(mag) / mag) * sk + ((1 - cos(mag)) / mag^2) * (sk * sk);
    else
        dC = eye(3) + sk + 0.5 * (sk * sk);
    end

    % Update DCM and normalize to maintain orthogonality: C <- 1.5*C - 0.5*C*(C'*C)
    Cs_e_next = Cs_e * dC;
    Cs_e_next = 1.5 * Cs_e_next - 0.5 * Cs_e_next * (Cs_e_next' * Cs_e_next);
end

% Normal gravity vector in ECEF frame: ge = C_g^e * gamma^g (Appendix A, F)
function ge = normal_gravity_ecef(x)
    % WGS-84 ellipsoid constants
    a = 6378137.0; f = 1 / 298.257223563; b = a * (1 - f);
    e2 = 2 * f - f^2; ep2 = (a^2 - b^2) / b^2;
    gamma_a = 9.7803253359; gamma_b = 9.8321849378;
    m = (7292115e-11)^2 * a^2 * b / 3986004.418e8;

    % Geodetic latitude, longitude, height (Bowring)
    p = hypot(x(1), x(2));
    th = atan2(x(3) * a, p * b);
    phi = atan2(x(3) + ep2 * b * sin(th)^3, p - e2 * a * cos(th)^3);
    lam = atan2(x(2), x(1));
    h = p / cos(phi) - a / sqrt(1 - e2 * sin(phi)^2);

    % Normal gravity on ellipsoid and at height h (Eq. 67, 69)
    s = sin(phi); c = cos(phi);
    g0 = (a * gamma_a * c^2 + b * gamma_b * s^2) / sqrt(a^2 * c^2 + b^2 * s^2);
    gh = g0 * (1 - (2 / a) * (1 + f + m - 2 * f * s^2) * h + (3 / a^2) * h^2);

    % North-south component (Eq. 70)
    f2 = -f + 2.5 * m + 0.5 * f^2 - (26 / 7) * f * m + 3.75 * m^2;
    gN = -((f2 - 0.5 * f^2 + 2.5 * f * m) / 6371000) * h * sin(2 * phi);

    % Transform from NED to ECEF: C_g^e * [gN; 0; gh] via Ce_g' (Appendix A, Eq. 48)
    ge = [-s * cos(lam) * gN - c * cos(lam) * gh;
          -s * sin(lam) * gN - c * sin(lam) * gh;
           c * gN            - s * gh           ];
end

% Direction Cosine Matrix from ECEF to local NED frame: Ce_g (Appendix A, Eq. 48)
function Ce_g = ecef2ned(x)
    a = 6378137.0; f = 1 / 298.257223563; b = a * (1 - f);
    ep2 = (a^2 - b^2) / b^2; e2 = 2 * f - f^2;
    p = hypot(x(1), x(2));
    th = atan2(x(3) * a, p * b);
    phi = atan2(x(3) + ep2 * b * sin(th)^3, p - e2 * a * cos(th)^3);
    lam = atan2(x(2), x(1));
    s = sin(phi); c = cos(phi);
    Ce_g = [-s * cos(lam), -s * sin(lam),  c;
            -sin(lam),      cos(lam),      0;
            -c * cos(lam), -c * sin(lam), -s];
end

% Skew-symmetric cross product matrix: skew(x) * y = cross(x, y)
function crossm = skew(x)
    crossm = [ 0,    -x(3),  x(2);
               x(3),  0,    -x(1);
              -x(2),  x(1),  0   ];
end

% Read IMU observations from binary file starting from epoch t0
function [dw, dv, t] = read_imu(filename, N, t0)
    fid = fopen(filename, 'rb');
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

% Compare estimated trajectory against reference trajectory in traj.txt
function compare_reference(t_est, x_est, v_est, traj_file)
    if ~exist(traj_file, 'file'), return; end
    traj = load(traj_file);
    t_ref = traj(:, 1);
    x_ref = traj(:, 2:4);
    v_ref = traj(:, 5:7);

    fprintf('\n--- Comparison against reference solution (%s) ---\n', traj_file);
    fprintf('%-11s | %-24s | %-24s\n', 'Time (s)', 'Pos Error [dX, dY, dZ] (m)', 'Vel Error [dVx, dVy, dVz] (m/s)');
    fprintf('%s\n', repmat('-', 1, 66));

    for i = 1:length(t_ref)
        [min_dt, idx] = min(abs(t_est - t_ref(i)));
        if min_dt < 0.003
            pos_err = x_est(idx, :) - x_ref(i, :);
            vel_err = v_est(idx, :) - v_ref(i, :);
            fprintf('%.4f | [%7.4f, %7.4f, %7.4f] | [%7.4f, %7.4f, %7.4f]\n', ...
                t_ref(i), pos_err, vel_err);
        end
    end
end