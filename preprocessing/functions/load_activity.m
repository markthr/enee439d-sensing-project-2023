function activity_data_struct = load_activity(subject_data_struct, activity)
    fn = subject_data_struct.sensors;
    activity_data_struct = struct;
    for i = 1:numel(fn)
        sensor = subject_data_struct.(fn{i});
        activity_data_struct.(fn{i}) = sensor.activity_data.(activity);
    end
    activity_data_struct.SubjectID = subject_data_struct.(fn{1}).subject;
    activity_data_struct.Activity = activity;
    activity_data_struct.sensors = fn;
end