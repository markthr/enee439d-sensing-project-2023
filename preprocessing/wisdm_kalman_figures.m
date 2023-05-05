%% Plot Interpolated Time Data for 29.A
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'A');
ds = align_sensor_times(ds, time_scale);


fig = figure;
fig.Position = fig_pos;
sensor = ds.p_gyr;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
Y = resample(X, t, fs, 'pchip', Dimension = 1);
components = ["X", "Y", "Z"];
N = 1000;
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(t(1:N), X(1:N, i), "k")
    scatter(t(1) + (0:(N-1))/fs, Y(1:N, i), "b.")
    xlim([t(1), t(1) + N/fs])
    ylabel('rad/s')
    subtitle(components(i) + " Component of " + sensor_names(4))
    hold off;
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')

%% Plot Orientation data based on Interpolation
fs = 20;
subject_data = load_subject(sensor_paths, 19);
ds = load_activity(subject_data, 'A');
ds = align_sensor_times(ds, time_scale);
% align time and make same length
X_acc_nu = xyz_to_mat(ds.w_acc);
X_gyr_nu = xyz_to_mat(ds.w_gyr);
t_acc_nu = double(ds.w_acc.TimeStampNanos)*1E-9;
t_gyr_nu = double(ds.w_gyr.TimeStampNanos)*1E-9;
X_acc = resample(X_acc_nu(2:(end-1), :), t_acc_nu(2:(end-1)), fs, 'pchip', Dimension = 1);

X_gyr = resample(X_gyr_nu, t_gyr_nu, fs, 'pchip', Dimension = 1);
fuse = imufilter('SampleRate',20);
%%
q = fuse(X_acc,X_gyr);
time = (0:size(X_acc,1)-1)/20;
plot(time,eulerd(q,'Z YX','frame'))
title('Orientation Estimate')
legend('Z-axis', 'Y-axis', 'X-axis')
xlabel('Time (s)')
ylabel('Rotation (degrees)')
%% Interpolated Time Delta for oversampled
subject_data = load_subject(sensor_paths, 1);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);
t = double(ds.w_acc.TimeStampNanos)*1E-9;
X = xyz_to_mat(ds.w_acc);
range = 401:900;
fig = figure;
fig.Position = fig_pos;
hold on;
plot(t(range), X(range, 3))

[y, ty] = resample(X(range, 3), t(range), fs, 'pchip');
plot(ty, y);
hold off;
%% Interpolated Time Delta for undersampled
subject_data = load_subject(sensor_paths, 50);
ds = load_activity(subject_data, 'M');
ds = align_sensor_times(ds, time_scale);
t = double(ds.w_acc.TimeStampNanos)*1E-9;
X = xyz_to_mat(ds.w_acc);
range = 2001:2900;
fig = figure;
fig.Position = fig_pos;
hold on;
plot(t(range), X(range, 3))

[y, ty] = resample(X(range, 3), t(range), fs, 'pchip');
plot(ty, y);
hold off;

%% functions
fs = 20;
thresh_diff = fs * 2 / 20;
% outputs the last index before a gap that exceeds the threshold
segment_ends = find(diff(t) > thresh_diff);

total_len = 0;
pred_len = 0;
i_start = 1; % first segment starts with the first index
segment_times = cell(length(segment_ends),1);
segment_values = cell(length(segment_ends),1);
for index = 1:numel(segment_ends)
    i_end = segment_ends(index);
    [Y, t_Y] = resample(X(i_start:i_end, :), t(i_start:i_end), fs, 'pchip');
    
    segment_times{index} = t_Y;
    segment_values{index} = Y;

    i_start = i_end;
end

% uniform sampled signal with gaps
X_ug = cat(1, segment_values{:});
t_ug = cat(1, segment_times{:});
%%
autoreg_len = 100;
autoreg_order = 40;
plot(t_ug, X_ug)
y = fillgaps(X_ug, autoreg_len, autoreg_order);
plot(y);

%%
function X = segmented_resampling(X, segments)

end
