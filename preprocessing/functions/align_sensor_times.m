function activity_data_struct = align_sensor_times(activity_data_struct, time_scale)
    % heuristic: assume remove any offset greater than 1 second in inital
    % times
    % heuristic: choose whatever ordering minimizes the average delay 
    % between beginning samples
    fn = activity_data_struct.sensors;
    n = numel(fn);

    % save t_0 of each sensor for determing which to place first
    beginnings = zeros(1,4, 'int64');
    % remove the part of offset that is greater than 1 second
    for i = 1:n
        activity_data_struct.(fn{i}).TimeStampNanos = ...
            remove_offset_from_nanos(activity_data_struct.(fn{i}).TimeStampNanos, time_scale);
        beginnings(i) = activity_data_struct.(fn{i}).TimeStampNanos(1);
    end
    
    % calculate resulting delays for all choices of t_0
    delays = zeros(n,n, 'int64');
    for i = 1:n
        delays(i,:) = beginnings - beginnings(i);
    end
    % add 1 second to negative values to get delay value
    mask = delays<0;
    delays(mask) =  time_scale + delays(mask);
    % select minimum average delay
    [~, index] = min(mean(delays,2));
    % generate offset to align the data
    alignments = zeros(1, n, 'int64') + beginnings - delays(index, :);

    % update time for each sensor
    for i = 1:n
        activity_data_struct.(fn{i}).TimeStampNanos = ...
            activity_data_struct.(fn{i}).TimeStampNanos - alignments(i);
    end
end

% simple helper function that removes any component of t_0 greater than 1
% second from a series of times in nanoseconds
function new_nanos = remove_offset_from_nanos(nanos, scale)
    new_nanos = nanos - idivide(nanos(1), scale)*scale;
end

