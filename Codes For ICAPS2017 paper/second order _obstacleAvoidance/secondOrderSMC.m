function     [posagents, velagents, Ck,Fox_old, Foy_old] = secondOrderSMC(posagents, velagents, Ck, muk, time, dt, DomainBounds, AgentForce, c, obstacles,Fox_old, Foy_old)
%% By: Hadi Salman
% 11/06/2016
global flag2 p limit
Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;
xmin = DomainBounds.xmin;
ymin = DomainBounds.ymin;

Nkx = size(muk, 1);
Nky = size(muk, 2);
Nagents = size(posagents, 1);

%% Updating Fourier Coefficients of Coverage Distribution
for iagent=1:Nagents
    xrel = posagents(iagent, 1) - xmin;
    yrel = posagents(iagent, 2) - ymin;
    for kx = 0:Nkx-1
        for ky = 0:Nky-1

%           hk = sqrt(Lx*Ly); % hadi
            hk = Lx*Ly; %hadi
            if kx ~= 0
                hk = hk * 0.5;
            end
            if ky ~= 0
                hk = hk * 0.5;
            end
            hk = sqrt(hk);%hadi

            Ck(kx+1, ky+1) = Ck(kx+1, ky+1) + (1/hk)*cos(kx * pi * xrel/Lx) * cos(ky * pi * yrel/Ly) * dt;
        end
    end
end

%% Computing controls
s=1.5;


%Ergodicity_Metric=0;
for iagent = 1:Nagents
    
    Bjx = 0.0;
    Bjy = 0.0;
    xrel = posagents(iagent, 1) - xmin;
    yrel = posagents(iagent, 2) - ymin;
    for kx = 0:Nkx-1
        for ky = 0:Nky-1
            lambda_k = 1.0 / ((1.0 + kx * kx + ky * ky)^s);
%             hk = sqrt(Lx*Ly); %hadi
              hk = Lx*Ly; %hadi
            if kx ~= 0
                hk = hk * 0.5;
            end
            if ky ~= 0
                hk = hk * 0.5;
            end
            hk = sqrt(hk); %hadi

            Bjx = Bjx + (lambda_k / hk) * (Ck(kx+1, ky+1) - Nagents*time*muk(kx+1, ky+1)) * (-kx * pi/Lx) * sin(kx * pi * xrel/Lx) * cos(ky * pi * yrel/Ly);
            Bjy = Bjy + (lambda_k / hk) * (Ck(kx+1, ky+1) - Nagents*time*muk(kx+1, ky+1)) * (-ky * pi/Ly) * cos(kx * pi * xrel/Lx) * sin(ky * pi * yrel/Ly);
        end
    end
    
    Bjnorm = sqrt(Bjx*Bjx + Bjy*Bjy);
    
    
%% OBSTACLE AVOIDANCE - Reference: Motion planning and collision avoidance using navigation vector fields
Fox = 0;
Foy = 0;
flag = 0;

if (flag2(iagent) == 0)%not changing as soon as it becomes in the vicinity of an obstacle
      p(iagent,:) =  ([Bjx,Bjy]./Bjnorm);
%     p(iagent,:) = - velagents(iagent, :)/(norm( velagents(iagent, :),2))
%       p(iagent,:) =  (c*velagents (iagent,:) + [Bjx, Bjy])/norm(c*velagents (iagent,:) + [Bjx, Bjy]);
end

flag2(iagent) = 0;
for iobstacle = 1:obstacles.number
ro = obstacles.r(iobstacle);%radius of obstacle
re = 0.01;%safety margin
r_robot= 0.00;% radius of robot

% obstacles' positions
xo = [obstacles.p(1,iobstacle)];
yo = [obstacles.p(2,iobstacle)];

 

% refer to the paper sec III
%valur of B at some point r > rz = ro+re+r
safety = 0.05;
r = ro + r_robot + re + safety; 
BF = ro^2 - r^2;
BZ = ro^2 - (ro + r_robot +re)^2;

       dx = posagents(iagent, 1)-xo;
       dy = posagents(iagent, 2)-yo;
       
       B = ro^2 - dx^2 - dy^2;
       
       if(B <= BF || B > 0)
          sigma = 1;% in vicinity of an obstacles
          flag2(iagent) = flag2(iagent) | 0;
       elseif (B > BF && B < BZ)
           flag2(iagent) = 1;
           flag = 1;%in vicinity of obstacle
%            a = 2/(BZ- BF)^3;
%            b = -3*(BZ+BF)/(BZ-BF)^3;
%            c = 6*BZ*BF/(BZ-BF)^3;
%            d = BZ^2*(BZ - 3*BF)/(BZ-BF)^3;
%            sigma =a*B^3 + b*B^2 + c*B + d;       
        sigma = 1 - abs((B - BF)/(BF-BZ));
       else 
           flag2(iagent) = 1;
           flag = 1;%in vicinity of obstacle
          sigma = 0;
       end
     if(flag == 1)  
       if (p(iagent,:)*[dx; dy] >= 0)
%             Fox =  (1-sigma)*(p(2)*(dx)*(dy) - p(1)*(dy)^2)/norm([p(2)*(dx)*(dy) - p(1)*(dy)^2,p(1)*(dx)*(dy) - p(2)*(dx)^2]) - sigma*((c*velagents(iagent,1) + Bjx)./norm(c*velagents (iagent,1) + Bjx));
%             Foy =  (1-sigma)*(p(1)*(dx)*(dy) - p(2)*(dx)^2)/norm([p(2)*(dx)*(dy) - p(1)*(dy)^2,p(1)*(dx)*(dy) - p(2)*(dx)^2]) - sigma*((c*velagents(iagent,2) + Bjy)./norm(c*velagents (iagent,2) + Bjy));
            Fox = Fox + (1-sigma)*(p(iagent,2)*(dx)*(dy) - p(iagent,1)*(dy)^2)/norm([p(iagent,2)*(dx)*(dy) - p(iagent,1)*(dy)^2,p(iagent,1)*(dx)*(dy) - p(2)*(dx)^2]) ;
            Foy = Foy + (1-sigma)*(p(iagent,1)*(dx)*(dy) - p(iagent,2)*(dx)^2)/norm([p(iagent,2)*(dx)*(dy) - p(iagent,1)*(dy)^2,p(iagent,1)*(dx)*(dy) - p(iagent,2)*(dx)^2]) ;
 
       end

       if (p(iagent,:)*[dx; dy]<0)
%             Fox =  (1-sigma)*(-p(iagent,1)*(dx).^2 - p(iagent,1)*(dy).^2)/norm([-p(iagent,1)*(dx).^2 - p(iagent,1)*(dy).^2,-p(iagent,2)*(dy).^2 - p(iagent,2)*(dx).^2])- sigma*((c*velagents(iagent,1) + Bjx)./norm(c*velagents (iagent,1) + Bjx));
%             Foy =  (1-sigma)*(-p(iagent,2)*(dy).^2 - p(iagent,2)*(dx).^2)/norm([-p(iagent,1)*(dx).^2 - p(iagent,1)*(dy).^2,-p(iagent,2)*(dy).^2 - p(iagent,2)*(dx).^2])- sigma*((c*velagents(iagent,2) + Bjy)./norm(c*velagents (iagent,2) + Bjy));
            Fox = Fox + (1-sigma)*(-p(iagent,1)*(dx).^2 - p(iagent,1)*(dy).^2)/norm([-p(iagent,1)*(dx).^2 - p(iagent,1)*(dy).^2,-p(iagent,2)*(dy).^2 - p(iagent,2)*(dx).^2]);
            Foy = Foy +(1-sigma)*(-p(iagent,2)*(dy).^2 - p(iagent,2)*(dx).^2)/norm([-p(iagent,1)*(dx).^2 - p(iagent,1)*(dy).^2,-p(iagent,2)*(dy).^2 - p(iagent,2)*(dx).^2]);
       end
     end
end
      if flag==1
       Fxdot = (Fox - Fox_old(iagent))/dt;
       Fydot = (Foy - Foy_old(iagent))/dt;
       Fox_old(iagent) = Fox;
       Foy_old(iagent) = Foy;
      
       L = sqrt(Fox^2 + Foy^2);
      end
       % Updating agent location based on SMC feedback control law
    % posagents(iagent,:) = posagents(iagent,:) + AgentSpeed * [Fox/L,Foy/L] * dt;
k=30;

%% Updating agent location based on SMC feedback control law
    posagents(iagent,:) = posagents(iagent, :) + velagents (iagent,:) * dt;
    if(flag == 1)
    velagents (iagent,:) = (velagents (iagent,:) + (k*(0.35 * [Fox/L,Foy/L] - velagents (iagent,:)) + .35*[Fxdot,Fydot])* dt);
    limit = [limit norm((k*(.35 * [Fox/L,Foy/L] - velagents (iagent,:)) + .35*[Fxdot,Fydot]))];
    else
    velagents (iagent,:) = (velagents (iagent,:) - AgentForce * ((c*velagents (iagent,:) + [Bjx, Bjy])./norm(c*velagents (iagent,:) + [Bjx, Bjy]))* dt);
    end
    
% Force the agent inward toward the domain with maximum force whenever they
% leave the domain
    velagents(iagent, :) = forcingInward(posagents(iagent,:), velagents(iagent, :), DomainBounds, AgentForce, dt);
end

end
%%
function [velagent] = forcingInward(posagent, velagent, DomainBounds, AgentForce, dt)

    xmin = DomainBounds.xmin;
    xmax = DomainBounds.xmax;
    ymin = DomainBounds.ymin;
    ymax = DomainBounds.ymax;

    if posagent(1) < xmin
        velagent(1) = velagent(1) + AgentForce*dt;
    end
    if posagent(1) > xmax
        velagent(1) = velagent(1) - AgentForce*dt;
    end
    if posagent(2) < ymin
        velagent(2) = velagent(2) + AgentForce*dt;    
    end
    if posagent(2) > ymax
        velagent(2) = velagent(2) - AgentForce*dt;    
    end

end
