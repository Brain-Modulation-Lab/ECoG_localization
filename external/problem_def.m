roll = 0;
pitch = 0; 
yaw = 0;
tx = -1.5;
ty = -2.5;
tz = -3;
alpha = 2;

Ro = angles2rotmat([pitch roll yaw]); 
RoSca = alpha*Ro; 
RoScaTra = [RoSca     [tx; ty; tz]; 
            0 0 0     1];