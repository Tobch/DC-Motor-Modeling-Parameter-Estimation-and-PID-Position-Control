% prepare_for_parameter_estimator.m
T = readtable('motor_seq_212336.csv');

% time vector (s)
t = T.t;

% supply voltage (use Vs column if present, otherwise assume 12V)
if any(strcmpi(T.Properties.VariableNames,'Vs'))
    Vs = T.Vs;
else
    Vs = 12*ones(size(t));
end

% duty -> fraction (0..1)
d = T.duty;
if max(d) <= 1
    duty = d;
elseif max(d) <= 100
    duty = d/100;
else
    duty = d/255;
end

% input (Volts)
u = Vs .* duty;

% output: convert RPM -> rad/s
y_rpm = T.omega;
y = y_rpm * 2*pi/60;

% sample time
Ts = median(diff(t));

% initial parameter guesses from your tfest output:
% Numerator = [0 170.586], Denominator = [1 20.1111]
num_tf = [0               0      73048.2923];
den_tf = [1        437.393      8614.5404];

% create iddata (recommended for Parameter Estimator)
z = iddata(y, u, Ts);
% give it a descriptive name in workspace
assignin('base','motor_data', z);

% also push initial vars to base workspace for the Simulink block
assignin('base','num_tf', num_tf);
assignin('base','den_tf', den_tf);

fprintf('Prepared workspace: variables motor_data (iddata), num_tf, den_tf, Ts\n');
