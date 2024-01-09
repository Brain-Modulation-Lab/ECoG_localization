function rigidmat = rotmat2rigidmat(rotmat, translation)
%UNTITLED Turn a 3x3 rotation matrix into a rigid transformation 4x4 matrix
if isempty(translation); translation = zeros(3, 1);

rigidmat = [rotmat    translation; 
            0 0 0     1];

end