function out = motor_log_and_ident(port, duration, varargin)

% ---------- parse inputs ----------
if nargin < 1 || strlength(string(port)) == 0, port = "COM3"; end
if nargin < 2 || isempty(duration),          duration = 20;   end

p = inputParser;
addParameter(p, "CSV", [], @(x) isempty(x) || (isscalar(x) && x>=1 && x<=1000));
addParameter(p, "File", "", @(s) isstring(s) || ischar(s));
addParameter(p, "SkipLog", false, @(b) islogical(b) || isnumeric(b));
addParameter(p, "UseFullSequence", false, @(b) islogical(b) || isnumeric(b));
parse(p, varargin{:});
csvMs          = p.Results.CSV;
outfile        = string(p.Results.File);
skipLog        = logical(p.Results.SkipLog);
useFullSeq     = logical(p.Results.UseFullSequence);

if outfile == ""
    outfile = "motor_seq_" + string(datestr(now,'HHMMSS')) + ".csv";
end

out = struct('T', table(), 'filename', "", 'params', struct(), ...
             'sys', [], 'quick_estimates', [NaN NaN NaN]);

% ---------- stage A: LOG (optional) ----------
if ~skipLog
    try
        logged = i_log_motor_sequence(port, duration, csvMs, outfile);
        out.T        = logged.T;
        out.filename = logged.filename;
    catch err
        fprintf(2, "✖ Logging failed: %s\n", err.message);
        return
    end
else
    % Load existing CSV
    assert(isfile(outfile), "CSV file '%s' not found.", outfile);
    T = readtable(outfile);
    out.T = T;
    out.filename = outfile;
end

% ---------- stage B: IDENTIFICATION (SI: rad/s per V) ----------
try
    ident = i_estimate_motor_params(out.filename, useFullSeq);
    out.params          = ident.params;
    out.sys             = ident.sys;
    out.quick_estimates = [ident.params.K_quick_rad_per_V, ...
                           ident.params.K_quick_RPM_per_V, ...
                           ident.params.tau_quick_s];
catch err
    fprintf(2, "✖ Identification failed: %s\n", err.message);
end

end

% ======================================================================
% ======================= INTERNAL SUBFUNCTIONS ========================
% ======================================================================

function out = i_log_motor_sequence(port, duration, csvMs, outfile)
baud = 115200;
s = [];
cleanupObj = onCleanup(@() i_safeClose(s)); %#ok<NASGU>

% ----- open serial -----
s = serialport(port, baud, "Timeout", 2);
configureTerminator(s, "LF");
flush(s);
pause(2.5);  % allow auto-reset

% drain boot chatter ~0.5 s
tDrain = tic;
while toc(tDrain) < 0.5 && s.NumBytesAvailable > 0
    readline(s);
end

% optional: set CSV interval on Arduino
if ~isempty(csvMs)
    writeline(s, "CSV " + string(csvMs));
    pause(0.05);
    flush(s);
end

% start sequence (confirm)
started = i_sendAndConfirm(s, "SEQ", "# sequence start", 2.0, 3);
if ~started
    warning("SEQ not acknowledged; trying D60 spin for 3 s, then SEQ once more...");
    writeline(s, "D60"); pause(3.0); writeline(s, "D0");
    started = i_sendAndConfirm(s, "SEQ", "# sequence start", 2.0, 1);
end
disp("✅ Logging...");

% capture <duration> seconds AFTER first valid row
cap  = max(ceil(duration*200), 200);
data = nan(cap,5); n = 0;

% wait for first numeric row (<=8 s)
firstFound = false; tStart = tic;
while ~firstFound && toc(tStart) < 8
    if s.NumBytesAvailable > 0
        line = strtrim(readline(s));
        if i_isGarbage(line), continue; end
        v = i_tryParse5(line);
        if ~isempty(v), n = 1; data(n,:) = v.'; firstFound = true; end
    else
        pause(0.002);
    end
end
if ~firstFound
    error("No valid numeric row within 8 s. Check Arduino is printing and motor has power.");
end

% timed capture
tStart = tic;
while toc(tStart) < duration
    if s.NumBytesAvailable > 0
        line = strtrim(readline(s));
        if i_isGarbage(line), continue; end
        v = i_tryParse5(line);
        if ~isempty(v)
            n = n + 1;
            if n > size(data,1)
                data = [data; nan(size(data,1),5)]; %#ok<AGROW>
            end
            data(n,:) = v.';
        end
    else
        pause(0.002);
    end
end

% stop + short drain
writeline(s, "STOP");
tDrain = tic;
while toc(tDrain) < 0.2 && s.NumBytesAvailable > 0, readline(s); end

% save + plot
data = data(1:n,:);
if isempty(data)
    warning("No data captured.");
    out = struct('T', table(), 'filename', "");
    return
end
T = array2table(data, 'VariableNames', {'t','duty','Vs','omega','theta'});

% --- NEW: compute u (input voltage) from duty and Vs ---
d = T.duty;
if max(d) <= 1.0
    duty = d;
elseif max(d) <= 100.0
    duty = d/100.0;
else
    duty = d/255.0;
end
T.u = T.Vs .* duty;   % volts

writetable(T, outfile);
fprintf("✅ Data saved: %s  (rows: %d)\n", outfile, height(T));

figure('Name','Motor Log');
subplot(3,1,1); plot(T.t, T.duty); ylabel('duty (%)'); grid on; title('Input');
subplot(3,1,2); plot(T.t, T.omega); ylabel('\omega (RPM)'); grid on; title('Speed');
subplot(3,1,3); plot(T.t, T.theta); ylabel('\theta (rev)'); xlabel('Time (s)'); grid on; title('Position');

out = struct('T', T, 'filename', outfile);
end

function out = i_estimate_motor_params(csvName, useFullSequence)
fprintf('>>> Loading %s\n', csvName);

T = readtable(csvName);
T.Properties.VariableNames = lower(T.Properties.VariableNames);
assert(all(ismember({'t','duty','omega'}, T.Properties.VariableNames)), ...
    'CSV must contain t, duty, omega columns.');

t = T.t(:);

% speed RPM -> rad/s
RPM2RAD = 2*pi/60;
y_rads  = T.omega(:) * RPM2RAD;

% optional position revolutions -> rad
theta_rad = [];
if ismember('theta', T.Properties.VariableNames)
    theta_rad = T.theta(:) * 2*pi;
end

% build input voltage from duty and Vs
if ismember('u', T.Properties.VariableNames)
    u = T.u(:);  % use already computed voltage
else
    if ismember('vs', T.Properties.VariableNames)
        Vsup = T.vs(:);
    else
        warning('Column "vs" missing; assuming Vs = 12 V.');
        Vsup = 12*ones(size(t));
    end
    d = T.duty(:);
    if max(d) <= 1.0
        duty = d;
    elseif max(d) <= 100.0
        duty = d/100.0;
    else
        duty = d/255.0;
    end
    u = Vsup .* duty;
end

Ts = median(diff(t));

% quick estimates (first step plateau) in SI
du   = diff(u);
idx0 = find(abs(du) > 0.5, 1, 'first'); % first step index
assert(~isempty(idx0), 'No clear input step detected.');

% next step change
idx1 = find(abs(du) > 0.5 & (1:numel(du))' > idx0, 1, 'first');
if isempty(idx1), idx1 = numel(t); end

plateau_start = idx0 + 1;
plateau_end   = idx1;
win_ss_start  = round(plateau_start + 0.8*(plateau_end - plateau_start));
win_ss_start  = max(plateau_start, min(win_ss_start, plateau_end-1));

y0  = mean(y_rads(max(1,idx0-20):idx0));
u0  = mean(u(max(1,idx0-20):idx0));
yss = mean(y_rads(win_ss_start:plateau_end-1));
uss = mean(u(win_ss_start:plateau_end-1));

K_si  = (yss - y0) / max(1e-9, (uss - u0));   % rad/s per V
K_rpm = K_si * (60/(2*pi));                   % RPM/V (for reference)

% tau via 10–90% exponential fit
rise_idx = (idx0+1):(idx1-1);
yspan = yss - y0;
tau  = NaN;
if ~isempty(rise_idx) && yspan > 0
    y10 = y0 + 0.10*yspan;
    y90 = y0 + 0.90*yspan;
    i10 = find(y_rads(rise_idx) >= y10, 1, 'first');
    i90 = find(y_rads(rise_idx) >= y90, 1, 'first');
    if ~isempty(i10) && ~isempty(i90) && i90 > i10
        tt = t(rise_idx(i10:i90)) - t(rise_idx(i10));
        yy = y_rads(rise_idx(i10:i90));
        z  = 1 - (yy - y0)/max(1e-9,yspan);
        z(z<=0) = eps;
        th = [ones(numel(tt),1) -tt] \ log(z);
        tau = 1/th(2);
    end
end

fprintf('Quick estimates (SI): K = %.3f rad/s/V  [= %.3f RPM/V], tau = %.4f s\n', ...
    K_si, K_rpm, tau);

% Identification (System Identification Toolbox: tfest/compare)
if useFullSequence
    fprintf('Using FULL sequence for identification (SI)...\n');
    zdata = iddata(y_rads, u, Ts);
else
    pre  = max(1, idx0 - round(0.5/Ts));
    post = min(length(t), idx0 + round(15/Ts));
    zdata = iddata(y_rads(pre:post), u(pre:post), Ts);
end

opt  = tfestOptions('InitialCondition','estimate','EnforceStability',true);
sys_1p0z     = tfest(zdata, 1, 0, opt);
try
    sys_1p0z_dly = tfest(zdata, 1, 0, NaN, opt);
catch
    sys_1p0z_dly = sys_1p0z;
end
sys_2p0z     = tfest(zdata, 2, 0, opt);

figure('Name','Measured vs Identified (SI: rad/s)');
compare(zdata, sys_1p0z, sys_1p0z_dly, sys_2p0z);
legend('Measured (rad/s)','1p0z','1p0z+delay','2p0z','Location','best');

fits = zeros(1,3);
[~, f] = compare(zdata, sys_1p0z);     fits(1) = f;
[~, f] = compare(zdata, sys_1p0z_dly); fits(2) = f;
[~, f] = compare(zdata, sys_2p0z);     fits(3) = f;

[bestFit, which] = max(fits);
models = {sys_1p0z, sys_1p0z_dly, sys_2p0z};
sys = models{which};
fprintf('Selected model: #%d  (Fit = %.2f%%)\n', which, bestFit);

[num, den] = tfdata(sys, 'v');
if isprop(sys,'IODelay'), L = sys.IODelay; else, L = 0; end
K_dc_id = dcgain(sys);   % rad/s per V

params = struct( ...
    'units','SI', ...
    'y_unit','rad/s', ...
    'theta_unit','rad', ...
    'K_quick_rad_per_V',K_si, ...
    'K_quick_RPM_per_V',K_rpm, ...
    'tau_quick_s',tau, ...
    'num',num,'den',den,'delay',L, ...
    'K_idc_rad_per_V',K_dc_id, ...
    'Ts',Ts,'csv',csvName);

if ~isempty(theta_rad), params.theta_rad = theta_rad; end

save('identified_params.mat','params');

fprintf('\nSimulink Transfer Fcn (OUTPUT = rad/s):\n');
fprintf('  Numerator   = [%s]\n', num2str(num));
fprintf('  Denominator = [%s]\n', num2str(den));
fprintf('  Transport Delay = %.4f s\n', L);
fprintf('  DC gain (identified) = %.3f rad/s per V  [= %.3f RPM/V]\n', ...
    K_dc_id, K_dc_id*(60/(2*pi)));

% Full-length comparison
figure('Name','Measured vs Identified (full, SI)');
zfull = iddata(y_rads, u, Ts);
compare(zfull, sys);
legend('Measured (rad/s)','Identified','Location','best');

out = struct('params', params, 'sys', sys);
end

% ---------- helpers ----------
function ok = i_sendAndConfirm(s, cmd, token, timeout_s, retries)
ok = false;
for k = 1:retries
    writeline(s, cmd);
    t0 = tic;
    while toc(t0) < timeout_s
        if s.NumBytesAvailable > 0
            line = strtrim(readline(s));
            if startsWith(line, token)
                ok = true; return;
            end
            v = i_tryParse5(line);
            if ~isempty(v), ok = true; return; end
        else
            pause(0.01);
        end
    end
end
end

function tf = i_isGarbage(line)
tf = (strlength(line)==0) || startsWith(line,"#") || ...
     startsWith(line,"t,duty,Vs,omega,theta",'IgnoreCase',true);
end

function v = i_tryParse5(line)
v = [];
vals = textscan(line,'%f%f%f%f%f','Delimiter',',','CollectOutput',true);
if ~isempty(vals) && ~isempty(vals{1}) && numel(vals{1})==5 && all(isfinite(vals{1}))
    v = vals{1};
else
    nums = sscanf(line,'%f,%f,%f,%f,%f');
    if numel(nums)==5 && all(isfinite(nums)), v = nums; end
end
end

function i_safeClose(s)
try
    if ~isempty(s)
        try, writeline(s,"STOP"); pause(0.05); catch, end
        try, flush(s); catch, end
    end
catch
end
end
