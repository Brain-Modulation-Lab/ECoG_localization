function T = cam2rigidtransmat(cam_pos, cam_target, cam_upvec)
b1 = normalize(cam_target(:) - cam_pos(:), "norm"); 


if ~isempty(cam_upvec)
    b3 = cam_upvec(:);
else
    % find the upvec closest to z direction
    z = [0 0 1]';
    b3 = z - dot(z,b1)*b1;
    b3 = normalize(b3, "norm");
end

b2 = cross(b1, -b3); 

assert(abs(dot(b1, b2)) < 10e-5); 
assert(abs(dot(b1, b3)) < 10e-5); 
assert(abs(dot(b2, b3)) < 10e-5); 

% T = [b1 b2 b3], "norm"); 
T = [[b1 b2 b3]         cam_target(:) ; 
     0 0 0               1          ];

% % testing 
% figure; 
% orig = repmat(cam_target', [3 1]); 
% quiver3(orig(:, 1), orig(:, 2), orig(:, 3), ...
%         T(1, 1:3)',  T(2, 1:3)',  T(3, 1:3)'); 


end