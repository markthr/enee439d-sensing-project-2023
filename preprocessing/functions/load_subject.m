function subject_data_struct = load_subject(sensor_paths_struct, subject_id)
    fn = fieldnames(sensor_paths_struct);
    subject_data_struct = struct;
    for i = 1:numel(fn)
        sensor_path = sensor_paths_struct.(fn{i});
        file_struct = sensor_path(subject_id);
        subject_data_struct.(fn{i}) = load([file_struct.folder '\' file_struct.name]);
    end

    subject_data_struct.SubjectID = subject_data_struct.(fn{1}).subject;
    subject_data_struct.sensors = fn;
end