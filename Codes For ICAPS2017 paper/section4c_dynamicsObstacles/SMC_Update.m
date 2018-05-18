%%BY: Hadi Salman
% OCT/09/2016
function [posagents, Ck] = SMC_Update(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed, obstacles)

Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;
xmin = DomainBounds.xmin;
ymin = DomainBounds.ymin;

Nkx = size(muk, 1);
Nky = size(muk, 2);
Nagents = size(posagents, 1);

% Updating Fourier Coefficients of Coverage Distribution
for iagent=1:Nagents
    xrel = posagents(iagent, 1) - xmin;
    yrel = posagents(iagent, 2) - ymin;
    for kx = 0:Nkx-1
        for ky = 0:Nky-1

             hk = Lx*Ly; 
            if kx ~= 0
                hk = hk * 0.5;
            end
            if ky ~= 0
                hk = hk * 0.5;
            end
            hk = sqrt(hk);

            Ck(kx+1, ky+1) = Ck(kx+1, ky+1) + (1/hk)*cos(kx * pi * xrel/Lx) * cos(ky * pi * yrel/Ly) * dt;
        end
    end
end

% Computing controls
s=1.5;

for iagent = 1:Nagents
Bjx = 0.0;
Bjy = 0.0;

    xrel = posagents(iagent, 1) - xmin;
    yrel = posagents(iagent, 2) - ymin;
    for kx = 0:Nkx-1
        for ky = 0:Nky-1
            lambda_k = 1.0 / ((1.0 + kx * kx + ky * ky)^s);

              hk = Lx*Ly; 
            if kx ~= 0
                hk = hk * 0.5;
            end
            if ky ~= 0
                hk = hk * 0.5;
            end
            hk = sqrt(hk); 

            Bjx = Bjx + (lambda_k / hk) * (Ck(kx+1, ky+1) - Nagents*time*muk(kx+1, ky+1)) * (-kx * pi/Lx) * sin(kx * pi * xrel/Lx) * cos(ky * pi * yrel/Ly);
            Bjy = Bjy + (lambda_k / hk) * (Ck(kx+1, ky+1) - Nagents*time*muk(kx+1, ky+1)) * (-ky * pi/Ly) * cos(kx * pi * xrel/Lx) * sin(ky * pi * yrel/Ly);
        end
    end

    Bjnorm = sqrt(Bjx*Bjx + Bjy*Bjy);

%% OBSTACLE AVOIDANCE - Reference: Motion planning and collision avoidance using navigation vector fields

Fox = 0;
Foy = 0;
product_sigma = 1;
%account for other agents
for other_agent=1:Nagents
   if other_agent ~= iagent
    ro = 0.01;%radius of obstacle
    re = 0.01;%safety margin
    r_robot= 0.00;% radius of robot

% other agents positions
    xo = [posagents(other_agent,1)];
    yo = [posagents(other_agent,2)];
    
  %direction of vector fields
    p = ([Bjx,Bjy]./Bjnorm);

% refer to the paper sec III
       smoothness = 0.01;
       r = ro + r_robot + re + smoothness; 
       BF = ro^2 - r^2;
       BZ = ro^2 - (ro + r_robot +re)^2;

       dx = posagents(iagent, 1)-xo;
       dy = posagents(iagent, 2)-yo;
       
       B = ro^2 - dx^2 - dy^2;
       
       if(B <= BF || B > 0)
          sigma = 1;
       elseif (B > BF && B < BZ)
%            a = 2/(BZ- BF)^3;
%            b = -3*(BZ+BF)/(BZ-BF)^3;
%            c = 6*BZ*BF/(BZ-BF)^3;
%            d = BZ^2*(BZ - 3*BF)/(BZ-BF)^3;
%            sigma =a*B^3 + b*B^2 + c*B + d;
           sigma =1-abs((B - BF)/(BF-BZ));
       else 
          sigma = 0;
       end
       product_sigma = product_sigma * sigma;
       if (p*[dx; dy] >= 0)
            Fox = Fox + (1-sigma)*(p(2)*(dx)*(dy) - p(1)*(dy)^2)/norm([p(2)*(dx)*(dy) - p(1)*(dy)^2,p(1)*(dx)*(dy) - p(2)*(dx)^2]);
            Foy = Foy + (1-sigma)*(p(1)*(dx)*(dy) - p(2)*(dx)^2)/norm([p(2)*(dx)*(dy) - p(1)*(dy)^2,p(1)*(dx)*(dy) - p(2)*(dx)^2]);
       end
       
       if (p*[dx; dy]<0)
            Fox =  Fox + (1-sigma)*(-p(1)*(dx).^2 - p(1)*(dy).^2)/norm([-p(1)*(dx).^2 - p(1)*(dy).^2,-p(2)*(dy).^2 - p(2)*(dx).^2]);
            Foy =  Foy +(1-sigma)*(-p(2)*(dy).^2 - p(2)*(dx).^2)/norm([-p(1)*(dx).^2 - p(1)*(dy).^2,-p(2)*(dy).^2 - p(2)*(dx).^2]);
       end
   end
end

% account for obstacles
for iobstacle = 1:obstacles.number

    ro = obstacles.r(iobstacle);%radius of obstacle
    re = 0.05;%safety margin
    r_robot= 0.00;% radius of robot

% obstacles' positions
    xo = [obstacles.p(1,iobstacle)];
    yo = [obstacles.p(2,iobstacle)];
    
  %direction of vector fields
    p = ([Bjx,Bjy]./Bjnorm);

% refer to the paper sec III
       smoothness = 0.1;
       r = ro + r_robot + re + smoothness; 
       BF = ro^2 - r^2;
       BZ = ro^2 - (ro + r_robot +re)^2;

       dx = posagents(iagent, 1)-xo;
       dy = posagents(iagent, 2)-yo;
       
       B = ro^2 - dx^2 - dy^2;
       
       if(B <= BF || B > 0)
          sigma = 1;
       elseif (B > BF && B < BZ)
%            a = 2/(BZ- BF)^3;
%            b = -3*(BZ+BF)/(BZ-BF)^3;
%            c = 6*BZ*BF/(BZ-BF)^3;
%            d = BZ^2*(BZ - 3*BF)/(BZ-BF)^3;
%            sigma =a*B^3 + b*B^2 + c*B + d;
           sigma =1-abs((B - BF)/(BF-BZ));
       else 
          sigma = 0;
       end
       product_sigma = product_sigma * sigma;
       if (p*[dx; dy] >= 0)
            Fox = Fox + (1-sigma)*(p(2)*(dx)*(dy) - p(1)*(dy)^2)/norm([p(2)*(dx)*(dy) - p(1)*(dy)^2,p(1)*(dx)*(dy) - p(2)*(dx)^2]);
            Foy = Foy + (1-sigma)*(p(1)*(dx)*(dy) - p(2)*(dx)^2)/norm([p(2)*(dx)*(dy) - p(1)*(dy)^2,p(1)*(dx)*(dy) - p(2)*(dx)^2]);
       end
       
       if (p*[dx; dy]<0)
            Fox =  Fox + (1-sigma)*(-p(1)*(dx).^2 - p(1)*(dy).^2)/norm([-p(1)*(dx).^2 - p(1)*(dy).^2,-p(2)*(dy).^2 - p(2)*(dx).^2]);
            Foy =  Foy +(1-sigma)*(-p(2)*(dy).^2 - p(2)*(dx).^2)/norm([-p(1)*(dx).^2 - p(1)*(dy).^2,-p(2)*(dy).^2 - p(2)*(dx).^2]);
       end
end
Fox = Fox - product_sigma * Bjx/Bjnorm;
Foy = Foy - product_sigma * Bjy/Bjnorm;
L = sqrt(Fox^2 + Foy^2);

    % Updating agent location based on SMC feedback control law
    posagents(iagent,:) = posagents(iagent,:) + AgentSpeed * [Fox/L,Foy/L] * dt;

    % reflecting agent in case it goes out of domain bounds
    posagents(iagent, :) = reflect_agent(posagents(iagent, :),DomainBounds);%hadi no Need for this

end

end

function [agent] = reflect_agent(agent, DomainBounds)

    xmin = DomainBounds.xmin;
    xmax = DomainBounds.xmax;
    ymin = DomainBounds.ymin;
    ymax = DomainBounds.ymax;

    if agent(1) < xmin
        agent(1) = xmin + (xmin - agent(1, 1));
    end
    if agent(1) > xmax
        agent(1) = xmax - (agent(1, 1) - xmax);
    end
    if agent(2) < ymin
        agent(2) = ymin + (ymin - agent(1, 2));
    end
    if agent(2) > ymax
        agent(2) = ymax - (agent(1, 2) - ymax);
    end

end
