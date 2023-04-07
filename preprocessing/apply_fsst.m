path = "J:/enee439d/datasets/wisdm-dataset";

watch_accel = dir(path + "/mat/watch/accel/*.mat");

%%
data_struct = load([watch_accel(1).folder '\' watch_accel(1).name]);
%%
data_a = data.activity_data.A;
t = double(data_a.TimeStampNanos) * 1E-9;
ts = duration(0,0,mean(t));
fsst(data_a.X,ts,'yaxis')
%%
spectrogram(data_a.X, 1000, 1000-1, 'yaxis')

%%
fs = 1000;
t = 0:1/fs:2;
ts = duration(0,0,1/fs);

x = chirp(t,100,1,200,'quadratic');

fsst(x,ts,'yaxis')

title('Quadratic Chirp')