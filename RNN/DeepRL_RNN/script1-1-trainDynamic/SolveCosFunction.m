function d = SolveCosFunction(params) % From Fan et al., 2024, Nature Human Bheavior
% params: [theta, scale]

% x_axis global
% y_axis local
% scale  local:global
% theta  angle between local and global
% [GHI  2er:         diff mid easy
%  DEF  1er:         diff mid easy
%  ABC] 0er/reward:  diff mid easy
theta = params(1);
scale = params(2);
%A, B, C
Ax=0; Ay=0;
Bx=1; By=0;
Cx=2; Cy=0;
%D, E, F
Dx = scale*cos(theta);
Dy = scale*sin(theta);
Ex = Dx+1;
Ey = Dy;
Fx = Dx+2;
Fy = Dy;

Gx = Dx*2;
Gy = Dy*2;
Hx = Gx+1;
Hy = Gy;
Ix = Gx+2;
Iy = Gy;
% X  = [Ax,Ay; Bx,By; Cx,Cy; Dx,Dy; Ex,Ey; Fx,Fy; Gx,Gy; Hx,Hy; Ix,Iy];
X  = [Ax,Ay; Dx,Dy; Gx,Gy; Bx,By; Ex,Ey; Hx,Hy; Cx,Cy; Fx,Fy; Ix,Iy];
d       = pdist(X,'euclidean');

end


