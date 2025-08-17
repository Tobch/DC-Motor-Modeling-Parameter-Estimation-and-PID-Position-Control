%% ===== PID Gain Calculation from Identified TF =====
% Load your identified transfer function parameters (from motor_log_and_ident)
% Example from your output:
num = [-1642.7752885709035 286714.952422509 1.3872860710096983E+6];      % Numerator of TF
den = [332.32922389216264 5824.1926660898635 17198.266704465146]; % Denominator of TF
K = num(end);                 % DC gain from numerator
a1 = den(2);
a0 = den(3);

% ===== Desired closed-loop poles =====
% Choose desired damping ratio and natural frequency
zeta = 1.2181;           % Overdamped
tr_desired = 0.05;  % desired rise time in seconds
wn = 2.2 / tr_desired; % wn ≈ 2.2 / tr for 0–100% rise

b1 = 2*zeta*wn;
b0 = wn^2;

% ===== PID calculation via coefficient matching =====
% Characteristic eq: s^2 + a1 s + a0 + K*(Kd s^2 + Kp s + Ki) = s^2 + b1 s + b0
Kd = (1 - 1) / K;       % simple approach: assume no derivative, can tune manually
Kp = (b1 - a1) / K;
Ki = (b0 - a0) / K;

% ===== Display results =====
fprintf('Calculated PID gains:\n');
fprintf('Kp = %.4f\n', Kp);
fprintf('Ki = %.4f\n', Ki);
fprintf('Kd = %.4f\n', Kd);

% Optional: create PID object for Simulink
PID_controller = pid(Kp, Ki, Kd);
