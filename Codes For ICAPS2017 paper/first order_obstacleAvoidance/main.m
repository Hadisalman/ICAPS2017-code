%% BY: Hadi Salman
%  10/09/2016
%% This codes prevents agent-agent collisions.(dynamic-obstacle avoidance)
%%JUST CLICK RUN! MAKE SURE YOU DON'T START IN AN OBSTACLE REGION
%% Setting domain bounds
DomainBounds.xmin = 0.0;
DomainBounds.xmax = 1.0;
DomainBounds.ymin = 0.0;
DomainBounds.ymax = 1.0;
Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;

%obstacles
obstacles.r = [0.1,0.05,0.08,0.05];% radii of obstacles
obstacles.p = [[0.5;0.5],[0.2;0.3],[0.7;0.2],[0.3;0.8]];%positions of obstacles[x1 y1
%                                                                               x2 y2
%                                                                               .....]
obstacles.number = numel(obstacles.r);

% Number of wave-numbers to be used
Nk = 50;%%

%% Calculating muk
res = 100;% resolution of discretization
xdel=Lx/res;
ydel=Ly/res;

mu=ones(res,res);
for xRange=0:xdel:Lx-xdel
         for yRange=0:ydel:Ly-ydel
            for i=1:obstacles.number
                if (xRange-obstacles.p(1,i))^2 + (yRange-obstacles.p(2,i))^2 <= obstacles.r(i)^2
                    mu(uint8(xRange*res+1),uint8(yRange*res+1)) = 0;
                end
            end
         end
end

mu=mu./sum(sum(mu));
[X,Y]=meshgrid(1:res,1:res);

muk=zeros(Nk,Nk);
Nkx = size(muk, 1);
Nky = size(muk, 2);
for kx = 0:Nkx-1
    for ky = 0:Nky-1
        
        hk=Lx*Ly; %using lim x->0 sinx/x=1
        if kx ~= 0
            hk = hk * 0.5;
        end
        if ky ~= 0
            hk = hk * 0.5;
        end
        hk = sqrt(hk);
        
        for xRange=0:xdel:Lx-xdel
            for yRange=0:ydel:Ly-ydel
                muk(kx+1, ky+1) = muk(kx+1, ky+1)+ mu(uint8(xRange*res+1),uint8(yRange*res+1)) *(1/hk)*cos(kx * pi * xRange/Lx) * cos(ky * pi * yRange/Ly);
                
            end
        end
        
        
    end
end

%% Initializing agent locations
Nagents = 4;
posagents = [0.35,0.7;0.05,0.3;0.1,0.95;0.83,0.03];%make sure not to start in an obstacle region!
AgentSpeed = 5;

colors = ['m','g','b','c'];

Nsteps = 5000;
dt = 0.001;

% Initializing Fourier coefficients of coverage distribution
Ck = zeros(Nk, Nk);

figure(1); hold on;
viscircles(obstacles.p',obstacles.r);
for xRange=0:xdel:Lx-xdel
         for yRange=0:ydel:Ly-ydel
            for i=1:obstacles.number
                if (xRange-obstacles.p(1,i))^2 + (yRange-obstacles.p(2,i))^2 <= obstacles.r(i)^2
                    scatter(xRange, yRange,2,'r','fill');
                end
            end
         end
end
axis equal;


ck_t = zeros(Nk, Nk);
phi_squared = zeros(Nsteps,1);
s= 1.5;


% Executing multiple steps of SMC algorithm
Ergodicity_Metric_save=0;
for it = 1:5000
    time = it * dt;
    [posagents, Ck] = SMC_Update(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed,obstacles);
    ck_t = Ck/(Nagents*time);
    for iagent = 1:Nagents
        plot(posagents(iagent, 1), posagents(iagent, 2), 'Color', colors(iagent) , 'Marker', 'o', 'MarkerSize', 1);
        axis([0 1 0 1])
        axis equal
        if it ~= 1
            h(iagent).Visible = 'off';
        end
        h(iagent) = scatter(posagents(iagent, 1), posagents(iagent, 2),colors(iagent),'fill');
%         pause(0.001)
    end
    
%     F(it) = getframe;
    it
    [Ergodicity_Metric] = Calculate_Ergodicity(Ck/Nagents/time, muk,DomainBounds);
    Ergodicity_Metric_save=[Ergodicity_Metric_save,Ergodicity_Metric];

end
%  F1=F(1:2500);
%  F2=F(2501:5000);
%  save('data.mat','F1','F2','Ergodicity_Metric_save','Nsteps');
%%
%  load('data.mat')
F=[F1,F2];
 movie(figure,F(5000),1,400);

%%
%  movie2avi(F(1:2500),'DynamicMultipleObstacles_ErgodicCoverage_NEW.avi','fps',100);
%%
close all
time=0:0.001:0.001*(Nsteps);
% figure;plot(time(2:end),Ergodicity_Metric_save(2:end))
figure;loglog(time(2:end),Ergodicity_Metric_save(2:end))
axis([0.001 5 0.0001,1])
xlabel('Time');
ylabel('Coverage Metric, \phi(t)');
% figure;loglog(time(2:end),Ergodicity_Metric_save(2:end))
% figure;plot(time(2:end),Ergodicity_Metric_save(2)./Ergodicity_Metric_save(2:end))
