function [posagents, Ck] = SMC_Update(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed)

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

%             hk = sqrt(Lx*Ly); %hadi
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

% Computing controls
s=1.5;
%Ergodicity_Metric=0;
for iagent = 1:Nagents
Bjx = 0.0; %hadi should be inside loop not outside
Bjy = 0.0; %hadi

    xrel = posagents(iagent, 1) - xmin;
    yrel = posagents(iagent, 2) - ymin;
    for kx = 0:Nkx-1
        for ky = 0:Nky-1
            lambda_k = 1.0 / ((1.0 + kx * kx + ky * ky)^s);
%               hk = sqrt(Lx*Ly); %hadi
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

    % Updating agent location based on SMC feedback control law
    posagents(iagent,:) = posagents(iagent,:) - AgentSpeed * ([Bjx,Bjy]./Bjnorm) * dt;
    
    %No Need for this
    % reflecting agent in case it goes out of domain bounds
    posagents(iagent, :) = reflect_agent(posagents(iagent, :), DomainBounds);
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
