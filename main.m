% INS GNSS Integration Framework
%
% Mechanization:
%   Coordinatized in the Earth Centered Earth Fixed (ECEF) frame (e-frame).
%   Attitude mechanization is strictly DCM-based (Direction Cosine Matrix, no quaternions).
%   Complementary / Error-State Extended Kalman Filter (15 states):
%     1:3   - ECEF position errors:     delta_x^e [m]
%     4:6   - ECEF velocity errors:     delta_v^e [m/s]
%     7:9   - ECEF attitude errors:     psi^e     [rad]
%     10:12 - Sensor accel bias errors: delta_b_f [m/s^2] (s-frame)
%     13:15 - Sensor gyro bias errors:  delta_b_w [rad/s] (s-frame)
%
% Generic notation:
% ECEF: Earth-Centered, Earth-Fixed frame (e-frame)
% NED:  Local level geodetic frame (North-East-Down, g-frame)
% Platform/Sensor: IMU body/sensor frame (s-frame = p-frame)
%
% Direction Cosine Matrices (DCM):
% Ce_g: ECEF to local level (NED) frame
% Cg_p: NED to platform frame (roll/pitch/yaw)
% Cs_p: Sensor to platform frame
% Cp_b: Platform to body frame (gyro stabilized IMU)
% Cp_c: Platform to camera frame (sensor orientation)
% Cm_c: Map to camera frame (omega/phi/kappa)

function main
    % System parameters and data filenames
    imu_filename  = 'data/primary_imu.dat';
    gnss_filename = 'data/ppp.txt';
    traj_filename = 'data/traj.txt';
    N = 50000;  % Number of IMU records to process

    % WGS-84 ellipsoid parameters
    global a b f e2 ep2 gamma_a gamma_b m
    a = 6378137.0;                                  % Semi-major axis [m]
    f = 1 / 298.257223563;                          % Flattening
    b = a * (1 - f);                                % Semi-minor axis [m]
    e2 = (a^2 - b^2) / a^2;                         % First eccentricity squared
    ep2 = (a^2 - b^2) / b^2;                        % Second eccentricity squared
    gamma_a = 9.7803253359;                         % Normal gravity at equator [m/s^2]
    gamma_b = 9.8321849378;                         % Normal gravity at pole [m/s^2]
    m = (7292115e-11)^2 * a^2 * b / 3986004.418e8;  % Angular velocity squared times gravitational parameter

    % =========================================================================
    % FIGURE 1 - BLOCK 1: Initial state
    % =========================================================================
    % Initial parameters from first line of data/traj.txt (epoch 198300.000259 s)
    t0   = 198300.000259;                              % Initial GPS Time of Week [s]
    x0_e = [3208319.5845; 633237.6578; 5458116.3716];  % Initial position in ECEF [m]
    v0_e = [-54.8084; -4.6903; 32.4978];               % Initial velocity in ECEF [m/s]
    Cg_p = [ 0.994637,  0.099168, -0.029369;           % NED to platform frame
            -0.102838,  0.918074, -0.382837;
            -0.011003,  0.383804,  0.923349];
    Cp_s = eye(3);                                      % Sensor to platform frame

    % Earth rotation rate vector and skew-symmetric matrix in ECEF [rad/s]
    wie_e = [0; 0; 7292115e-11];
    Wie_e = skew(wie_e);

    % Initial attitude DCM: Cs_e0 = (Cp_s * Cg_p * Ce_g0)' (Appendix A, Eq. 48)
    Ce_g0 = ecef2ned(x0_e);
    Cs_e0 = (Cp_s * Cg_p * Ce_g0)';
    Cs_e0 = normalize(Cs_e0);

    % Initial navigation states
    x    = x0_e;
    v    = v0_e;
    Cs_e = Cs_e0;

    % Initial sensor bias estimates (Eq. 47a-b)
    b_f = [0; 0; 0];  % Accelerometer bias in sensor frame [m/s^2]
    b_w = [0; 0; 0];  % Gyroscope bias in sensor frame [rad/s]

    % GNSS Antenna Lever-Arm in body/sensor frame [m] (from data/L390_20260331.mnt)
    % Primary antenna position relative to IMU reference point:
    ls = [-2.0210; 0.5093; -1.0411];

    % Tuning Parameters
    sigma_p0    = [0.1; 0.1; 0.1];                 % Initial position std [m]
    sigma_v0    = [0.05; 0.05; 0.05];              % Initial velocity std [m/s]
    sigma_psi0  = [0.5; 0.5; 1.0] * (pi / 180);    % Initial attitude roll/pitch/heading std [rad]
    sigma_ba0   = [0.01; 0.01; 0.01];              % Initial accel bias std [m/s^2] (~1000 ug)
    sigma_bg0   = [0.01; 0.01; 0.01] * (pi / 180); % Initial gyro bias std [rad/s] (~0.5 deg/s)

    P = diag([sigma_p0.^2; sigma_v0.^2; sigma_psi0.^2; sigma_ba0.^2; sigma_bg0.^2]);

    % Filter stochastic process noise tuning parameters (Power Spectral Densities):
    filter_params.vrw = 0.05;                      % Velocity Random Walk (accel white noise) [m/s / sqrt(s)]
    filter_params.arw = 0.005 * (pi / 180);        % Angular Random Walk (gyro white noise) [rad / sqrt(s)]
    filter_params.s_ba = 1e-4;                     % Accel bias random walk PSD [(m/s^2)^2 / s]
    filter_params.s_bg = 1e-6 * (pi / 180);        % Gyro bias random walk PSD [(rad/s)^2 / s]
    filter_params.tau_a = Inf;                     % Accel bias Gauss-Markov correlation time [s] (Inf = Random Walk)
    filter_params.tau_g = Inf;                     % Gyro bias Gauss-Markov correlation time [s] (Inf = Random Walk)

    % Load Sensor Observations
    % Read high-rate IMU observations (delta-v and delta-theta)
    fprintf('Reading IMU observations from %s...\n', imu_filename);
    [dw, dv, t_imu] = read_imu(imu_filename, N, t0);
    fprintf('IMU records: %d | Start ToW: %.6f s | End ToW: %.6f s\n', ...
        size(dw, 1), t_imu(1), t_imu(end));

    % Read aiding GNSS observations (PPP fixes from data/ppp.txt)
    fprintf('Reading GNSS aiding fixes from %s...\n', gnss_filename);
    gnss = read_gnss_ppp(gnss_filename, t_imu(1) - 1.0, t_imu(end) + 1.0);
    fprintf('Loaded %d GNSS fixes in the dataset time window.\n', length(gnss.t));

    % Pre-allocate trajectory history arrays
    x_all      = zeros(N, 3);
    v_all      = zeros(N, 3);
    b_f_all    = zeros(N, 3);
    b_w_all    = zeros(N, 3);
    P_diag_all = zeros(N, 15);

    x_all(1, :)      = x';
    v_all(1, :)      = v';
    b_f_all(1, :)    = b_f';
    b_w_all(1, :)    = b_w';
    P_diag_all(1, :) = diag(P)';

    % Index tracker for aiding observations
    next_gnss_idx = 1;
    t_last_kf = t_imu(1);
    num_updates = 0;

    fprintf('\n====================================================================\n');
    fprintf('Starting Aided INS Processing (Figure 1 Flowchart, DCM Mechanization)\n');
    fprintf('====================================================================\n');
    tic;

    for k = 1:N-1
        dt = t_imu(k+1) - t_imu(k);
        if dt <= 0, dt = 0.005; end                      % Time step [s] (IMU rate)
        t_curr = t_imu(k+1);

        % =====================================================================
        % FIGURE 1 - BLOCK 2: IMU measurements
        % =====================================================================
        % Raw sensor increments: dw (angular increment [rad]), dv (velocity increment [m/s])
        dw_raw = dw(k, :)';
        dv_raw = dv(k, :)';

        % Correct raw increments using currently estimated sensor biases:
        dw_k = dw_raw - b_w * dt;
        dv_k = dv_raw - b_f * dt;

        % =====================================================================
        % FIGURE 1 - BLOCK 3: Integrate navigation state
        % =====================================================================
        % 1. DCM attitude integration (Eqs. 31, 61) via matrix exponential
        [Cs_e_next, rho] = integrate_dcm(Cs_e, dw_k, dt, wie_e);

        % 2. Normalize DCM to preserve orthonormality
        Cs_e_next = normalize(Cs_e_next);

        % 3. Numerical integration of position and velocity (Appendix C.2, Eqs. 59-62b)
        [x, v, ~] = integrate_pos_vel(x, v, dv_k, rho, dt, Cs_e_next, Wie_e);
        Cs_e = Cs_e_next;

        % =====================================================================
        % FIGURE 1 - DECISION: Aiding obs. available?
        % =====================================================================
        % Check if a GNSS fix is available within tolerance around the current epoch
        aiding_available = false;
        if next_gnss_idx <= length(gnss.t)
            t_gnss = gnss.t(next_gnss_idx);
            % Tolerance is half an IMU epoch (dt / 2)
            if abs(t_curr - t_gnss) <= (dt / 2 + 1e-6)
                aiding_available = true;
                gnss_epoch_idx = next_gnss_idx;
                next_gnss_idx = next_gnss_idx + 1;
            elseif t_gnss < t_curr - (dt / 2)
                % Skip overdue GNSS fix
                next_gnss_idx = next_gnss_idx + 1;
            end
        end

        % If GNSS epochs do not coincide exactly with IMU samples, consider
        % state interpolation/extrapolation (Appendix G, page 14).

        if aiding_available
            % Time elapsed since last Kalman filter update
            dt_kf = t_curr - t_last_kf;
            if dt_kf <= 0, dt_kf = dt; end

            % Specific force in sensor frame [m/s^2]
            f_s = dv_k / dt;

            % =================================================================
            % FIGURE 1 - BLOCK 5: Time update error state cov. matrix
            % =================================================================
            [P, ~, ~] = time_update_cov(P, x, v, Cs_e, f_s, dt_kf, Wie_e, filter_params);

            % =================================================================
            % FIGURE 1 - BLOCK 6: Measurement update
            % =================================================================
            % GNSS Position fix and measurement covariance
            z_pos = gnss.pos(gnss_epoch_idx, :)';
            R_pos = gnss.cov_mat{gnss_epoch_idx};

            % 1. Compensate for lever arm: x_ant = x + Cs_e * ls
            % 2. Innovation: dz = z_pos - x_ant
            % 3. Design matrix: H = [I_3, 0_3x3, -skew(Cs_e * ls), 0_3x3, 0_3x3]
            % 4. Kalman Gain: K = P * H' / (H * P * H' + R_pos)
            % 5. Error state: dx_hat = K * dz
            % 6. Posterior covariance: P = (I - K*H) * P
            [dx_hat, P, ~, ~, ~] = measurement_update_gnss(x, Cs_e, P, z_pos, R_pos, ls);

            % =================================================================
            % FIGURE 1 - BLOCK 7: Correct navigation state
            % =================================================================
            % Correct position: x = x + dx_hat(1:3) (Eq. 46a)
            % Correct velocity: v = v + dx_hat(4:6) (Eq. 46b)
            % Correct attitude: Cs_e = (I + skew(psi)) * Cs_e (Eq. 46c, strictly DCM)
            % Re-normalize Cs_e
            % Correct sensor biases: b_f = b_f + dx_hat(10:12) (Eq. 47a)
            %                        b_w = b_w + dx_hat(13:15) (Eq. 47b)
            % Reset error state: dx_hat = 0 (closed-loop complementary filter)
            % -----------------------------------------------------------------
            [x, v, Cs_e, b_f, b_w] = correct_nav_state(x, v, Cs_e, b_f, b_w, dx_hat);

            t_last_kf = t_curr;
            num_updates = num_updates + 1;
        end

        % Store navigation estimates for logging and trajectory analysis
        x_all(k+1, :)      = x';
        v_all(k+1, :)      = v';
        b_f_all(k+1, :)    = b_f';
        b_w_all(k+1, :)    = b_w';
        P_diag_all(k+1, :) = diag(P)';
    end
    t_elapsed = toc;

    fprintf('Processing complete: %d epochs in %.3f s | GNSS updates performed: %d\n', ...
        N, t_elapsed, num_updates);
    fprintf('Final estimated position [m]:   [%.4f, %.4f, %.4f]\n', x_all(end, :));
    fprintf('Final estimated velocity [m/s]: [%.4f, %.4f, %.4f]\n', v_all(end, :));
    fprintf('Final accel bias [mg]:          [%.2f, %.2f, %.2f]\n', b_f_all(end, :) * (1e3 / 9.80665));
    fprintf('Final gyro bias [deg/h]:        [%.2f, %.2f, %.2f]\n', b_w_all(end, :) * (180 / pi * 3600));

    % =========================================================================
    % Validation & Benchmark against Reference Trajectory
    % =========================================================================
    compare_reference(t_imu, x_all, v_all, traj_filename);

    % Plot trajectories and diagnostics in local geodetic frame (NED)
    plot_trajectory_g(t_imu, x_all, x0_e, Ce_g0, traj_filename);
    plot_filter_diagnostics(t_imu, b_f_all, b_w_all, P_diag_all);
end


% =============================================================================
% LOCAL FUNCTIONS
% =============================================================================

% -----------------------------------------------------------------------------
% Time update error state cov. matrix
% -----------------------------------------------------------------------------
function [Pp, Phi, Q] = time_update_cov(P, x, v, Cs_e, f_s, dt, Wie_e, params)
    % Evaluates continuous system error dynamics matrix F (Eq. 45)
    F = error_system_matrix_F(x, v, Cs_e, f_s, Wie_e, params);

    % State transition matrix Phi
    % For aiding intervals (e.g. dt = 1.0 s), matrix exponential expm(F*dt)
    % (Eq. 5) or interval sub-stepping (App. D, Eq. 63) ensures numerical stability.
    % Truncated Taylor series (Eq. 6) may also be evaluated.
    Phi = expm(F * dt);

    % Process noise covariance Q
    % Construct continuous noise gain matrix G and spectral density W
    G = zeros(15, 12);
    G(4:6,   1:3)   = -Cs_e;        % Accel white noise into velocity error
    G(7:9,   4:6)   = -Cs_e;        % Gyro white noise into attitude error
    G(10:12, 7:9)   = eye(3);       % Accel bias driving noise
    G(13:15, 10:12) = eye(3);       % Gyro bias driving noise

    W = diag([repmat(params.vrw^2, 3, 1); ...
              repmat(params.arw^2, 3, 1); ...
              repmat(params.s_ba, 3, 1); ...
              repmat(params.s_bg, 3, 1)]);

    % 1st-order approximation: Q = G * W * G' * dt (Eq. 12 / 64)
    Q = G * W * G' * dt;

    % Prior covariance: P_prior = Phi * P * Phi' + Q
    Pp = Phi * P * Phi' + Q;

    % Enforce exact symmetry
    %Pp = 0.5 * (Pp + Pp');
end

% -----------------------------------------------------------------------------
% Continuous system error dynamic matrix F
% -----------------------------------------------------------------------------
function F = error_system_matrix_F(x, v, Cs_e, f_s, Wie_e, params)
    % Specific force in ECEF frame [m/s^2]
    fe = Cs_e * f_s;

    % Gravitational gradient tensor Gamma_e (Appendix E, Eqs. 65-66)
    Gamma_e = gravitational_gradient_tensor(x);

    % Construct 15x15 continuous F matrix according to Eq. (45)
    F = zeros(15, 15);

    % Row 1-3: Position error derivative: d(delta_x)/dt = delta_v (Eq. 37)
    F(1:3, 4:6) = eye(3);

    % Row 4-6: Velocity error derivative (Eq. 36)
    % d(delta_v)/dt = (Gamma_e - Wie_e*Wie_e)*delta_x - 2*Wie_e*delta_v - [fe x]*psi - Cs_e*delta_b_f
    F(4:6, 1:3)   = Gamma_e - Wie_e * Wie_e;
    F(4:6, 4:6)   = -2 * Wie_e;
    F(4:6, 7:9)   = -skew(fe);
    F(4:6, 10:12) = -Cs_e;   % Sensor accel bias into velocity error (Eq. 45)

    % Row 7-9: Attitude error derivative (Eq. 44)
    % d(psi)/dt = -Wie_e*psi - Cs_e*delta_b_w
    F(7:9, 7:9)   = -Wie_e;
    F(7:9, 13:15) = -Cs_e;   % Sensor gyro bias into attitude error (Eq. 45)

    % Row 10-12: Accelerometer bias dynamics (First-order Gauss-Markov or Random Walk)
    if isfinite(params.tau_a) && params.tau_a > 0
        F(10:12, 10:12) = -(1 / params.tau_a) * eye(3);
    end

    % Row 13-15: Gyroscope bias dynamics (First-order Gauss-Markov or Random Walk)
    if isfinite(params.tau_g) && params.tau_g > 0
        F(13:15, 13:15) = -(1 / params.tau_g) * eye(3);
    end
end

% -----------------------------------------------------------------------------
% Gravitational gradient tensor in ECEF frame
% -----------------------------------------------------------------------------
function Gamma_e = gravitational_gradient_tensor(x)
    % ECEF position unit vector e = x / ||x||
    r = norm(x);
    if r < 1e3, r = 6371000; end
    e = x / r;

    % Schuler frequency: omega_s = sqrt(||g|| / Re) (Eq. 66)
    Re = 6378137.0;       % Earth radius [m]
    g_mag = 9.81;         % Approximate gravity magnitude [m/s^2]
    omega_s = sqrt(g_mag / Re);

    % Gravitational gradient tensor: Gamma_e = omega_s^2 * (3*e*e' - I) (Eq. 65)
    Gamma_e = (omega_s^2) * (3 * (e * e') - eye(3));
end

% -----------------------------------------------------------------------------
% Measurement Update
% -----------------------------------------------------------------------------
function [dx_hat, P, dz, S, K] = measurement_update_gnss(x, Cs_e, Pp, z_pos, R_pos, ls)
    % Transform antenna lever-arm from sensor to ECEF frame: l^e = Cs_e * l^s
    le = Cs_e * ls;

    % Predicted observation: position of GNSS antenna in ECEF
    x_ant = x + le;

    % Measurement innovation
    dz = z_pos - x_ant;

    % Design matrix H
    H = zeros(3, 15);
    H(1:3, 1:3) = eye(3);
    H(1:3, 7:9) = -skew(le);

    % Innovation covariance: S = H * P * H' + R (Eq. 17d)
    S = H * Pp * H' + R_pos;

    % Kalman Gain K
    K = (Pp * H') / S;

    % Error state estimate
    dx_hat = K * dz;

    % Posterior covariance
    % Joseph stabilized form: P = (I - K*H)*P*(I - K*H)' + K*R*K'
    IKH = eye(15) - K * H;
    P = IKH * Pp * IKH' + K * R_pos * K';
    %P = 0.5 * (P + P');
end

% -----------------------------------------------------------------------------
% Correct navigation state 
% -----------------------------------------------------------------------------
function [x, v, Cs_e, b_f, b_w] = correct_nav_state(x, v, Cs_e, b_f, b_w, dx_hat)
    % Position correction
    x = x + dx_hat(1:3);

    % Velocity correction
    v = v + dx_hat(4:6);

    % Attitude DCM correction
    % Cs_e = (I + Psi_e) * Cs_e, where Psi_e = skew(psi)
    psi = dx_hat(7:9);
    Cs_e = (eye(3) + skew(psi)) * Cs_e;

    % Re-orthogonalize / normalize DCM to maintain mathematical orthogonality
    Cs_e = normalize(Cs_e);

    % Sensor bias corrections
    b_f = b_f + dx_hat(10:12);
    b_w = b_w + dx_hat(13:15);
end

% -----------------------------------------------------------------------------
% Numerical integration of position & velocity
% -----------------------------------------------------------------------------
function [x, v, dv_e] = integrate_pos_vel(x, v, dv, rho, dt, Cs_e_next, Wie_e)
    % Midpoint average DCM over update interval
    Cbar = Cs_e_next * (eye(3) - 0.5 * skew(rho));

    % Normal gravity vector in ECEF frame: ge = C_g^e * gamma^g (Appendix F)
    ge = normal_gravity_ecef(x);

    % Eq. (59): Velocity increment transformed to ECEF with Coriolis and gravity
    dv_e = Cbar * dv - (2 * Wie_e * v - ge) * dt;

    % Eq. (62a-b): Position and velocity updates
    x = x + v * dt + dv_e * (dt / 2);
    v = v + dv_e;
end

% -----------------------------------------------------------------------------
% DCM attitude integration
% -----------------------------------------------------------------------------
function [Cs_e_next, rho] = integrate_dcm(Cs_e, dw, dt, wie_e)
    % Rotation increment relative to ECEF frame
    rho = dw - Cs_e' * wie_e * dt;
    mag = norm(rho);
    Gamma = skew(rho);

    % Relative DCM increment (matrix exponential via Rodrigues formula)
    if mag > 1e-12
        dC = eye(3) + (sin(mag) / mag) * Gamma + ((1 - cos(mag)) / mag^2) * (Gamma * Gamma);
    else
        dC = eye(3) + Gamma + 0.5 * (Gamma * Gamma);
    end

    % Update DCM: Cs_e(t_{k+1}) = Cs_e(t_k) * dC (Eq. 31)
    Cs_e_next = Cs_e * dC;
end

% -----------------------------------------------------------------------------
% Normalize a DCM to be orthogonal: C = 1.5*C - 0.5*C*(C'*C)
% -----------------------------------------------------------------------------
function C = normalize(C)
    C = 1.5 * C - 0.5 * C * (C' * C);
end

% -----------------------------------------------------------------------------
% Read GNSS PPP observations from text file (e.g., data/ppp.txt)
% -----------------------------------------------------------------------------
function gnss = read_gnss_ppp(filename, t_start, t_end)
    if nargin < 2, t_start = -Inf; end
    if nargin < 3, t_end   = Inf;  end

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open GNSS file: %s', filename);
    end

    % Read columns: ToW, X, Y, Z, Vx, Vy, Vz, stdX, stdY, stdZ, cov21, cov31, cov32
    raw = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f', 'CommentStyle', '#');
    fclose(fid);

    all_t   = raw{1};
    all_pos = [raw{2}, raw{3}, raw{4}];
    all_vel = [raw{5}, raw{6}, raw{7}];
    all_std = [raw{8}, raw{9}, raw{10}];
    all_cov = [raw{11}, raw{12}, raw{13}];

    % Filter within requested time window
    valid_mask = (all_t >= t_start) & (all_t <= t_end);
    gnss.t   = all_t(valid_mask);
    gnss.pos = all_pos(valid_mask, :);
    gnss.vel = all_vel(valid_mask, :);
    gnss.std = all_std(valid_mask, :);
    gnss.cov = all_cov(valid_mask, :);

    % Build 3x3 covariance matrix cell array for fast lookup
    M = length(gnss.t);
    gnss.cov_mat = cell(M, 1);
    for i = 1:M
        sx2 = gnss.std(i, 1)^2;
        sy2 = gnss.std(i, 2)^2;
        sz2 = gnss.std(i, 3)^2;
        c21 = gnss.cov(i, 1);
        c31 = gnss.cov(i, 2);
        c32 = gnss.cov(i, 3);
        gnss.cov_mat{i} = [sx2, c21, c31; ...
                           c21, sy2, c32; ...
                           c31, c32, sz2];
    end
end

% -----------------------------------------------------------------------------
% ECEF coordinates to geodetic coordinates (latitude, longitude, height)
% -----------------------------------------------------------------------------
function [lat, lon, h] = ecef2geod(x)
    global a b f e2 ep2 gamma_a gamma_b m
    p = hypot(x(1), x(2));
    my = atan2(x(3) * a, p * b);
    lat = atan2(x(3) + ep2 * b * sin(my)^3, p - e2 * a * cos(my)^3);
    lon = atan2(x(2), x(1));
    h = p / cos(lat) - a / sqrt(1 - e2 * sin(lat)^2);
end

% -----------------------------------------------------------------------------
% Normal gravity vector in ECEF frame: ge = C_g^e * gamma^g
% -----------------------------------------------------------------------------
function ge = normal_gravity_ecef(x)
    global a b f e2 ep2 gamma_a gamma_b m
    [lat, ~, h] = ecef2geod(x);
    Ce_g = ecef2ned(x);

    g0 = (a * gamma_a * cos(lat)^2 + b * gamma_b * sin(lat)^2) / ...
         sqrt(a^2 * cos(lat)^2 + b^2 * sin(lat)^2);
    gh = g0 * (1 - (2 / a) * (1 + f + m - 2 * f * sin(lat)^2) * h + (3 / a^2) * h^2);

    f2 = -f + 2.5 * m + 0.5 * f^2 - (26 / 7) * f * m + 3.75 * m^2;
    gN = -((f2 - 0.5 * f^2 + 2.5 * f * m) / 6371000) * h * sin(2 * lat);

    ge = Ce_g' * [gN; 0; gh];
end

% -----------------------------------------------------------------------------
% Direction Cosine Matrix from ECEF to local NED frame
% -----------------------------------------------------------------------------
function Ce_g = ecef2ned(x)
    [lat, lon, ~] = ecef2geod(x);
    Ce_g = [-sin(lat) * cos(lon), -sin(lat) * sin(lon),  cos(lat);   % North
            -sin(lon),              cos(lon),               0;       % East
            -cos(lat) * cos(lon), -cos(lat) * sin(lon), -sin(lat)];  % Down
end

% -----------------------------------------------------------------------------
% Skew-symmetric cross product matrix: skew(x) * y = cross(x, y) 
% -----------------------------------------------------------------------------
function crossm = skew(x)
    crossm = [ 0,    -x(3),  x(2);
               x(3),  0,    -x(1);
              -x(2),  x(1),  0   ];
end

% -----------------------------------------------------------------------------
% Read IMU observations from binary file starting from epoch t0
% -----------------------------------------------------------------------------
function [dw, dv, t] = read_imu(filename, N, t0)
    fid = fopen(filename, 'rb');
    if fid == -1
        error('Could not open IMU file: %s', filename);
    end

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

    dw = dw / 206264806.247; % Gyro angular increments [rad]
    dv = dv / 20000000.0;     % Accel velocity increments [m/s]
end

% -----------------------------------------------------------------------------
% Compare estimated trajectory against reference solution in traj.txt
% -----------------------------------------------------------------------------
function compare_reference(t_est, x_est, v_est, traj_file)
    if ~exist(traj_file, 'file'), return; end
    traj = load(traj_file);
    t_ref = traj(:, 1);
    x_ref = traj(:, 2:4);
    v_ref = traj(:, 5:7);

    fprintf('\n--- Comparison against reference solution (%s) ---\n', traj_file);
    fprintf('%-11s | %-26s | %-26s\n', 'Time (s)', 'Pos Error [dX, dY, dZ] (m)', 'Vel Error [dVx, dVy, dVz] (m/s)');
    fprintf('%s\n', repmat('-', 1, 72));

    pos_errors = [];
    vel_errors = [];
    for i = 1:length(t_ref)
        [min_dt, idx] = min(abs(t_est - t_ref(i)));
        if min_dt < 0.003
            pe = x_est(idx, :) - x_ref(i, :);
            ve = v_est(idx, :) - v_ref(i, :);
            pos_errors = [pos_errors; pe];
            vel_errors = [vel_errors; ve];
            if mod(i, 10) == 1 || i == length(t_ref)
                fprintf('%.4f | [%7.3f, %7.3f, %7.3f] | [%7.3f, %7.3f, %7.3f]\n', ...
                    t_ref(i), pe, ve);
            end
        end
    end

    if ~isempty(pos_errors)
        rmse_pos = sqrt(mean(pos_errors.^2, 1));
        rmse_vel = sqrt(mean(vel_errors.^2, 1));
        fprintf('%s\n', repmat('-', 1, 72));
        fprintf('RMSE Position (X, Y, Z) [m]:   [%.3f, %.3f, %.3f] | Total 3D: %.3f m\n', ...
            rmse_pos, norm(rmse_pos));
        fprintf('RMSE Velocity (Vx, Vy, Vz) [m/s]: [%.3f, %.3f, %.3f] | Total 3D: %.3f m/s\n', ...
            rmse_vel, norm(rmse_vel));
    end
end

% -----------------------------------------------------------------------------
% Plot estimated position in local geodetic g-frame (NED)
% -----------------------------------------------------------------------------
function plot_trajectory_g(t, x_est, x0_e, Ce_g0, traj_file)
    if nargin < 5, traj_file = ''; end

    x_est_g = (x_est - x0_e') * Ce_g0';
    time = t - t(1);

    has_ref = false;
    if ~isempty(traj_file) && exist(traj_file, 'file')
        traj = load(traj_file);
        t_ref = traj(:, 1);
        idx_ref = t_ref >= t(1) - 0.1 & t_ref <= t(end) + 0.1;
        if any(idx_ref)
            x_ref = traj(idx_ref, 2:4);
            x_ref_g = (x_ref - x0_e') * Ce_g0';
            t_ref_sub = t_ref(idx_ref) - t(1);
            has_ref = true;
        end
    end

    % Figure 1: Horizontal Trajectory in NED
    figure('Name', 'Horizontal Trajectory in g-frame (NED)', 'NumberTitle', 'off');
    plot(x_est_g(:, 2), x_est_g(:, 1), 'b-', 'LineWidth', 1.5);
    hold on;
    if has_ref
        plot(x_ref_g(:, 2), x_ref_g(:, 1), 'k--', 'LineWidth', 1.2);
    end
    plot(x_est_g(1, 2), x_est_g(1, 1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(x_est_g(end, 2), x_est_g(end, 1), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    title('Aided INS Horizontal Trajectory in g-frame (NED)');
    if has_ref
        legend('Aided INS', 'Reference', 'Start', 'End', 'Location', 'northeast');
    else
        legend('Aided INS', 'Start', 'End', 'Location', 'northeast');
    end

    % Figure 2: Position Components in NED vs Time
    figure('Name', 'Position Components in g-frame (NED)', 'NumberTitle', 'off');
    labels = {'North [m]', 'East [m]', 'Down [m]'};
    for ax = 1:3
        subplot(3, 1, ax);
        plot(time, x_est_g(:, ax), 'b-', 'LineWidth', 1.5);
        hold on; grid on;
        if has_ref
            plot(t_ref_sub, x_ref_g(:, ax), 'k--', 'LineWidth', 1.2);
        end
        ylabel(labels{ax});
        if ax == 1
            title('Aided INS Position Components in g-frame (NED)');
            if has_ref, legend('Aided INS', 'Reference', 'Location', 'northeast'); end
        end
    end
    xlabel('Time [s]');
    drawnow;
end

% -----------------------------------------------------------------------------
% Plot filter diagnostics: Estimated sensor biases and 3-sigma uncertainty bounds
% -----------------------------------------------------------------------------
function plot_filter_diagnostics(t, b_f_est, b_w_est, P_diag_est)
    time = t - t(1);

    figure('Name', 'Estimated Sensor Biases & Uncertainties', 'NumberTitle', 'off');

    % Accelerometer biases [mg]
    subplot(2, 2, 1);
    plot(time, b_f_est * (1e3 / 9.80665), 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Accel Bias [mg]');
    title('Estimated Accelerometer Biases');
    legend('X (s-frame)', 'Y (s-frame)', 'Z (s-frame)', 'Location', 'northeast');

    % Gyroscope biases [deg/h]
    subplot(2, 2, 2);
    plot(time, b_w_est * (180 / pi * 3600), 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Gyro Bias [deg/h]');
    title('Estimated Gyroscope Biases');
    legend('X (s-frame)', 'Y (s-frame)', 'Z (s-frame)', 'Location', 'northeast');

    % Position 3-sigma uncertainty [m]
    subplot(2, 2, 3);
    sigma_pos = sqrt(max(0, P_diag_est(:, 1:3)));
    plot(time, 3 * sigma_pos, 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('3\sigma Position Bound [m]');
    title('Position Uncertainty Bounds (3\sigma)');
    legend('X_e', 'Y_e', 'Z_e', 'Location', 'northeast');

    % Velocity 3-sigma uncertainty [m/s]
    subplot(2, 2, 4);
    sigma_vel = sqrt(max(0, P_diag_est(:, 4:6)));
    plot(time, 3 * sigma_vel, 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('3\sigma Velocity Bound [m/s]');
    title('Velocity Uncertainty Bounds (3\sigma)');
    legend('V_{x,e}', 'V_{y,e}', 'V_{z,e}', 'Location', 'northeast');

    drawnow;
end