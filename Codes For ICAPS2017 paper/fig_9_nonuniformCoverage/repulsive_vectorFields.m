%%BY: Hadi Salman
% OCT/05/2016
%% Reference: Motion planning and collision avoidance using navigation vector fields
% JUST CLICK RUN!

%%

ro = 1;%radius of obstacle
re = 0.0;%safety margin
r_robot= 0.0;% radius of robot

figure;
hold on
viscircles([-2,-2],ro);
axis equal;
axis([-5 2 -5 2]);
% scatter(0,0,'x');
% text(0,0,'  target');

% obstacles' positions
xo = [-2];
yo = [-2];
N = numel(xo); % number of obstacles

% repulsive vector field (for one obstacle)
phi = atan2(- yo, -xo)+ pi;
p = [cos(phi), sin(phi)];% direction from obstacle to origin
% p = [1,0];
resolution = 0.3;
[x,y] = meshgrid(-5.0:resolution:2.0,-5.0:resolution:2.0);

% %repalsive vector fields around obstacle o
% %lambda = 1: p.deltaR > 0 (Check ref. paper sec: III)
% Fox1 = p(2)*(x-xo).*(y-yo) - p(1)*(y-yo).^2;
% Foy1 = p(1)*(x-xo).*(y-yo) - p(2)*(x-xo).^2;
% 
% %lambda = 0 p.deltaR < 0  (Check ref. paper sec: III)
% Fox0 = -p(1)*(x-xo).^2 - p(1)*(y-yo).^2;
% Foy0 = -p(2)*(y-yo).^2 - p(2)*(x-xo).^2;

%%
% refer to the paper sec III
%valur of B at some point r > rz = ro+re+r
r = ro + r_robot + re + 1; 
BF = ro^2 - r^2;
BZ = ro^2 - (ro + r_robot +re)^2;

i=0;
j=0;
Fox = zeros(size(x));
Foy = zeros(size(x));
for xi = -5.0:resolution:2.0
   j=j+1;
   i=0;
    for yi = -5.0:resolution:2.0
       i=i+1;
       dx = xi-xo;
       dy = yi-yo;
       
       B = ro^2 - dx^2 - dy^2;
       
       if(B <= BF || B > 0)
          sigma = 1;
       elseif (B > BF && B < BZ)
           a = 2/(BZ- BF)^3;
           b = -3*(BZ+BF)/(BZ-BF)^3;
           c = 6*BZ*BF/(BZ-BF)^3;
           d = BZ^2*(BZ - 3*BF)/(BZ-BF)^3;
           sigma =a*B^3 + b*B^2 + c*B + d;
       else 
          sigma = 0;
       end
           
       if (p*[dx; dy] >= 0)
            Fox(i, j) =  (1-sigma)*(p(2)*(dx)*(dy) - p(1)*(dy)^2)/norm([(p(2)*(dx)*(dy) - p(1)*(dy)^2),(p(1)*(dx)*(dy) - p(2)*(dx)^2)]);
            Foy(i, j) =  (1-sigma)*(p(1)*(dx)*(dy) - p(2)*(dx)^2)/norm([(p(2)*(dx)*(dy) - p(1)*(dy)^2),(p(1)*(dx)*(dy) - p(2)*(dx)^2)]);
       end
       
       if (p*[dx; dy]<0)
            Fox(i, j) =  (1-sigma)*(-p(1)*(dx).^2 - p(1)*(dy).^2)/norm([(-p(1)*(dx).^2 - p(1)*(dy).^2),(-p(2)*(dy).^2 - p(2)*(dx).^2)]);
            Foy(i, j) =  (1-sigma)*(-p(2)*(dy).^2 - p(2)*(dx).^2)/norm([(-p(1)*(dx).^2 - p(1)*(dy).^2),(-p(2)*(dy).^2 - p(2)*(dx).^2)]);
       end
    end
end

% Fox = Fox1;
% Foy =  Foy1;
% figure(1)
L = sqrt(Fox.^2 + Foy.^2);
quiver(x,y,Fox./L,Foy./L);

