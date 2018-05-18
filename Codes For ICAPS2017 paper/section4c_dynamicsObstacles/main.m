%% BY: Hadi Salman
%  10/09/2016
%%This codes prevents agent-agent collisions.(dynamic-obstacle avoidance)
%%JUST CLICK RUN! MAKE SURE YOU DON'T START IN AN OBSTACLE REGION
%% Setting domain bounds
DomainBounds.xmin = 0.0;
DomainBounds.xmax = 1.0;
DomainBounds.ymin = 0.0;
DomainBounds.ymax = 1.0;
Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;

% Number of wave-numbers to be used
Nk = 50;%%

%% Calculating muk
res = 100;% resolution of discretization
xdel=Lx/res;
ydel=Ly/res;

mu=ones(res,res);

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

%obstacles

obstacles.r = [0.1,0.05,0.05];% radii of obstacles
obstacles_pos = @(it)[[0.1*cos(2*pi*it/3000)+0.5;0.2*sin(2*pi*it/1000)+0.5], ...
                      [0.2;0.2*sin(2*pi*it/2000)+0.6],....
                      [0.2*sin(2*pi*it/1000)+0.3;0.2]];%,[0.2;0.3],[0.7;0.2],[0.3;0.8]positions of obstacles[x1 y1
%                                                                               x2 y2
%                                                                               .....]
obstacles.number = numel(obstacles.r);

figure(1);hold on

ck_t = zeros(Nk, Nk);
phi_squared = zeros(Nsteps,1);
s= 1.5;


% Executing multiple steps of SMC algorithm
Ergodicity_Metric_save=0;
for it = 1:Nsteps
    time = it * dt;
    obstacles.p = obstacles_pos(it);
    [posagents, Ck] = SMC_Update(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed,obstacles);
    ck_t = Ck/(Nagents*time);
    for iagent = 1:Nagents
%         plot(posagents(iagent, 1), posagents(iagent, 2), 'Color', colors(iagent) , 'Marker', 'o', 'MarkerSize', 1);
        axis([0 1 0 1])
        axis equal
        xlim([0,1])
        ylim([0,1])
        if it ~= 1
            h(iagent).Visible = 'off';
        end
        h(iagent) = scatter(posagents(iagent, 1), posagents(iagent, 2),colors(iagent),'fill');
        viscircles(obstacles_pos(it)',obstacles.r); 
        pause(0.001)
    end
    
    if mod(it,100)==0
       fprintf('Iteration: %i/%i  \n', it,Nsteps) 
    end
    [Ergodicity_Metric] = Calculate_Ergodicity(Ck/Nagents/time, muk,DomainBounds);
    Ergodicity_Metric_save=[Ergodicity_Metric_save,Ergodicity_Metric];

end

%% Plotting the metric of ergodicity
close all
time=0:0.001:0.001*(Nsteps);
figure;loglog(time(2:end),Ergodicity_Metric_save(2:end))
axis([0.001 5 0.0001,1])
xlabel('Time');
ylabel('Coverage Metric, \phi(t)');
