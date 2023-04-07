path = "J:/enee439d/datasets/wisdm-dataset";

watch_accel = dir(path + "/mat/watch/accel/*.mat");

%%
data_struct = load([watch_accel(1).folder '\' watch_accel(1).name]);
%% compare FSST with and without HPF
data_a = data_struct.activity_data.A;
t = double(data_a.TimeStampNanos) * 1E-9;
ts = duration(0,0,mean(t));
[s,f,t] = fsst(data_a.X,ts);
ax = subplot(2, 1, 1);
waterfall(ax, f,t,abs(s)'.^2)
% apply LPF to remove gravity
X_hp = highpass(data_a.X, 2E-4);
[s,f,t] = fsst(X_hp,ts);
ax = subplot(2, 1, 2);
waterfall(ax, f,t,abs(s)'.^2)
%% iterate over fields and apply FFT
transform_fields = ['X', 'Y', 'Z'];
activities = fieldnames(data_struct.activity_data);

for i = 1:numel(activities)
    activity = activities{i};
    disp(activity)
    activity_data = data_struct.activity_data.(activity);
    % calculate average sampling time
    ts_avg = mean(diff(activity_data.TimeStampNanos));
    for j = 1:numel(transform_fields)
        transform_field = transform_fields(j);
        disp(transform_field)
        
        raw = activity_data.(transform_field);
        [s, f, t] = stft(raw);
    end
    
    
end
%%
fs = 1000;
t = 0:1/fs:2;
ts = duration(0,0,1/fs);

x = chirp(t,100,1,200,'quadratic');

fsst(x,ts,'yaxis')

title('Quadratic Chirp')