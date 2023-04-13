%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";

window_time = 5; % how often to process a new window in seconds
max_window_error = 0.1; % the maximum percent a window can be too short by 
window_time_shift = 1;
nufft_length = 100;
% Setting window_memory = 0 indicates that only the current window should be included for processing.

% get letter for the 18 activities, A-S but without N
activities = char([1:13, 15:19] + 'A' - 1);

time_scale = int64(1E9); % seconds to nanos
fs = 20;
%% load directories
watch_accel = dir(path + "/mat/watch/accel/*.mat");
watch_gyro = dir(path + "/mat/watch/gyro/*.mat");
phone_accel = dir(path + "/mat/phone/accel/*.mat");
phone_gyro = dir(path + "/mat/phone/gyro/*.mat");

%% check data to ensure subjects line up
for i = 1:numel(watch_accel)
    get_subj = "^data_(\d{4})_[a-z]+_[a-z]+\.mat$";
    subject = regexp(watch_accel(i).name, get_subj, 'tokens');
    match_subj = "^data_" + subject{1}{1} + "_[a-z]+_[a-z]+\.mat$";
    assert(regexp(watch_gyro(i).name, match_subj) == 1, "Subject IDs do not match");
    assert(regexp(phone_accel(i).name, match_subj) == 1, "Subject IDs do not match");
    assert(regexp(phone_gyro(i).name, match_subj) == 1, "Subject IDs do not match");
end
%%
ds = load_data_struct(watch_accel(1));
%% compare FSST with and without HPF
data_a = ds.activity_data.A;
t = double(data_a.TimeStampNanos) * 1E-9;
ts = mean(diff(t));
fs = 1/ts;
X_hp = highpass(data_a.X, 2E-4);
%%
[s,f,t] = fsst(data_a.X,fs);
ax = subplot(2, 1, 1);
waterfall(ax, f,t,abs(s)'.^2)
% apply LPF to remove gravity

[s,f,t] = fsst(X_hp,fs);
ax = subplot(2, 1, 2);
waterfall(ax, f,t,abs(s)'.^2)

%% fsst vs stft
subplot(2, 1, 1);
stft(X_hp,fs, 'FrequencyRange', "onesided")
% apply LPF to remove gravity
subplot(2, 1, 2);
fsst(X_hp,fs, 'yaxis');

%% make plots using nustft
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    [s2, t] = nustft(X, double(ds.(fn{i}).TimeStampNanos)*1E-9, fs, window_time, window_time_shift, max_window_error);
    f = (0:49)/100*20;
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle(['Subject: ', int2str(ds.(fn{1}).SubjectID(1)), ', Activity: ', activities(activity_index)])
%% make table
tabulate(ds.p_acc.Group)
%% plot time delta
seconds_w_acc = double(ds.w_acc.TimeStampNanos(400:900))*1E-9;
plot(401:900, diff(seconds_w_acc))
ylim([0, 0.4])
xlabel('Sample')
ylabel('Time Difference')
%% a
Y = nufft(double(ds.p_acc.hh_vec), double(ds.p_acc.TimeStampNanos), [], 2);
disp(size(Y))
for group_index = 1:numel(ds.p_acc.Group(end))
    unfft(ds.p_gyr.hh_vec(i,:))
    
end

%%

%%
activity_data = watch_accel_data.activity_data.A;
ts_avg = mean(diff(activity_data.TimeStampNanos));
w = kaiser(ceil(10E9/ts_avg),10);
raw = activity_data.Y;
[s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
subplot(2,1,1)
surf(t, f, 20*log10(abs(s).^2), 'EdgeColor', 'none');
activity_data = watch_accel_data.activity_data.B;
ts_avg = mean(diff(activity_data.TimeStampNanos));
raw = activity_data.Y;
[s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
subplot(2,1,2)
surf(t, f, 20*log10(abs(s).^2), 'EdgeColor', 'none');
%%
fs = 1000;
t = 0:1/fs:2;
ts = duration(0,0,1/fs);

x = chirp(t,100,1,200,'quadratic');

fsst(x,ts,'yaxis')

title('Quadratic Chirp')
%%
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    [s2, t] = nustft(X, double(ds.(fn{i}).TimeStampNanos)*1E-9, fs, window_time, window_time_shift, max_window_error);
    f = (0:49)/100*20;
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')

end
sgtitle(['Subject: ', int2str(ds.(fn{1}).SubjectID(1)), ', Activity: ', activities(activity_index)])

%% function declarations
function data_struct = load_data_struct(file_struct)
    data_struct = load([file_struct.folder '\' file_struct.name]);
end

function new_nanos = remove_offset_from_nanos(nanos, scale)
    new_nanos = nanos - idivide(nanos(1), scale)*scale;
end


function alignments = get_alignment(beginnings, time_scale)
    % heuristic: choose whatever ordering minimizes the average delay 
    % between beginning samples
    n = length(beginnings);
    % calculate resulting delays for all choices of t_0
    delays = zeros(n,n, 'int64');
    for i = 1:n
        delays(i,:) = beginnings - beginnings(i);
    end
    % add 1 second to negative values to get delay value
    mask = delays<0;
    % select minimum average delay
    delays(mask) = time_scale - delays(mask);
    [~, index] = min(mean(delays,2));
    alignments = zeros(n,n, 'int64') - beginnings(index);
    alignments(mask(index, :)) = alignments(mask(index, :)) + time_scale;
    
end

function matrix = xyz_to_mat(struct, sequence)
    if(~exist('sequence','var'))
        matrix = [struct.X; struct.Y; struct.Z];
    else
        matrix = [struct.X(sequence); struct.Y(sequence); struct.Z(sequence)];
    end
    
end

function [starts, ends, targets] = partition_time_sequence(t, t_win, t_err, nufft_length)
    start_of_sec = zeros(uint32(floor(t(end))), 1, 'uint32');
    starts = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    ends = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    targets = 5:(numel(start_of_sec)-1);
    % get initial value for start times
    for i = 1:numel(start_of_sec)
        start_of_sec(i) = find(t>=(i-1), 1, 'first');
    end
    % fix sequences to be a multiple of 100 in length
    for i = (t_win+1):(numel(start_of_sec))
        i_start = start_of_sec(i-t_win);
        i_end = start_of_sec(i) - 1;
        seq_len = i_end - i_start + 1;
        mod_nufft_len = mod(seq_len, nufft_length);
        if(mod_nufft_len == 0)
            % the window lined up perfectly
        elseif(seq_len < nufft_length)
            % the window was too short, extend it
            if(i_start - (nufft_length-mod_nufft_len) <= 0)
                % cannot extend backwards, extend forwards. 
                % 
                % This complicates analysis, but this corner case can happen 
                % only once per activity so it should hopefully not be major.
                % TODO: if the data towards the end gets discard, might as well
                % shift everything forward such that there is enough data
                % before the first window that backwards extension works
                i_end = i_end + (nufft_length-mod_nufft_len);
            else
                % can extend backwards
                i_start = i_start - (nufft_length-mod_nufft_len);
            end
        % if execution reaches here, sequence must be too long
        elseif(t(i_end) - t(i_start+mod_nufft_len) > t_win - t_err)
            % Window can be shortened without making the time covered by the
            % window too short.
            i_start = i_start + mod_nufft_len;
        else
            % the window is too long, but must be further lengthened to cover
            % the minimum desired timespan while having the samples be a factor
            % of 100
            if(i_start - (nufft_length-mod_nufft_len) <= 0)
                % cannot extend backwards, extend forwards. Same drawbacks as
                % above hold.
                i_end = i_end + (nufft_length-mod_nufft_len);
            else
                % can extend backwards
                i_start = i_start - (nufft_length-mod_nufft_len);
            end
        end
        starts(i-t_win) = i_start;
        ends(i-t_win) = i_end;
    end
end

% TODO: current, only handles overlap times of multiples of 1 second
function [s2, t] = nustft(x, t, fs, window_time, window_time_shift, max_window_error)
    nufft_length = 100;
    [starts, ends, targets] = partition_time_sequence(t, window_time, window_time * max_window_error, nufft_length);
    disp(t(ends(1:3)))
    [x_dim, ~] = size(x);
    window_indices = 1:window_time_shift:numel(starts);
    s2 = zeros(nufft_length/2, x_dim, nufft_length/2);
    for i = window_indices
        n = ends(i) - starts(i) + 1;
        f = double(0:(n/2-1))/double(n)*fs;
        Y = nufft(x(:,starts(i):ends(i)), t(starts(i):ends(i)), f, 2);
        bin_size = idivide(ends(i)-starts(i)+1, nufft_length);
        if(bin_size == 1)
            s2(i,:,:) = abs(Y .* conj(Y));
        else
            s2(i,:,:) = squeeze(sum(reshape(abs(Y.*conj(Y)),bin_size, x_dim, [])))./double(bin_size);
        end
    end
    t = targets;
end
