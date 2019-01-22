function [x,y] = circle_tube(R,cx,cy)
%%%%%%%%%%%%%%%%%%%
% 画圆函数
%%%%%%%%%%%%%%%%%%%
nb_pts = 10;
alpha=0:pi/nb_pts:2*pi;%角度[0,2*pi]
%R=2;%半径
x=R*cos(alpha)+cx;
y=R*sin(alpha)+cy;
% plot(cx,cy,'g.',x,y,'g.');
plot( x,y,'k.');
grid on;
% axis equal;