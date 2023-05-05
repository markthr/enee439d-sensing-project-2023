function matrix = xyz_to_mat(struct, sequence)
    if(~exist('sequence','var'))
        matrix = [struct.X, struct.Y, struct.Z];
    else
        matrix = [struct.X(sequence), struct.Y(sequence), struct.Z(sequence)];
    end
end