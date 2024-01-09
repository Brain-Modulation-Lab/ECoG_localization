function T = cam2rigidtransmat(cam_pos, cam_target, cam_upvec)
b1 = cam_target(:) - cam_pos(:); 
b3 = cam_upvec(:);
b2 = cross(b1, -b3); 
assert(abs(dot(b1, b2)) < 10e-5); 
assert(abs(dot(b1, b3)) < 10e-5); 
assert(abs(dot(b2, b3)) < 10e-5); 

T = normalize([b1 b2 b3], "norm"); 
T = [T(1:3, 1:3)         cam_target ; 
     0 0 0               1          ];
end