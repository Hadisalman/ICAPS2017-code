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
                                                                               .....]
obstacles.number = numel(obstacles.r);

% Number of wave-numbers to be used
Nk = 50;%%
livePlot = true; %set <true> if you want live plot for trajectories, other <false> for faster execution
%% Calculating muk
res = 100;% resolution of discretization
mu=ones(res,res);
xdel=Lx/res;
ydel=Ly/res;

% for xRange=0:xdel:Lx-xdel
%          for yRange=0:ydel:Ly-ydel
%             for i=1:obstacles.number
%                 if (xRange-obstacles.p(1,i))^2 + (yRange-obstacles.p(2,i))^2 <= obstacles.r(i)^2
%                     mu(uint8(xRange*res+1),uint8(yRange*res+1)) = 0;
%                 end
%             end
%          end
% end
mu=mu./sum(sum(mu));
[X,Y]=meshgrid(1:res,1:res);
surf(X,Y,mu);
close;
% muk = dct2(mu,Nk,Nk);
%For Matlab DCT to work: (k1,k2)=1,2,..(L1-L2)

muk=zeros(Nk,Nk);
xdel=Lx/res;
ydel=Ly/res;
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
ck_t = zeros(Nk, Nk);
figure(1); hold on;
viscircles(obstacles.p',obstacles.r);
%color inside the circles
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

phi_squared = zeros(Nsteps,1);
s= 1.5;
% Executing multiple steps of SMC algorithm
Ergodicity_Metric_save=0;
tic
for it = 1:Nsteps
    time = (it) * dt;
    [posagents, Ck] = SMC_Update(posagents, Ck, muk, time, dt, DomainBounds, AgentSpeed);
    ck_t = Ck/(Nagents*time);
    for iagent = 1:Nagents
        plot(posagents(iagent, 1), posagents(iagent, 2), 'Color', colors(iagent), 'Marker', 'o', 'MarkerSize', 1);
        if livePlot == true
                pause(0.001)
        end
    end
    if mod(it,100)==0
       fprintf('Iteration: %i/%i  \n', it,Nsteps) 
    end
    
    [Ergodicity_Metric] = Calculate_Ergodicity(Ck/Nagents/time, muk,DomainBounds);
    Ergodicity_Metric_save=[Ergodicity_Metric_save,Ergodicity_Metric];
    
end
toc
%% Plotting the metric of ergodicity 
time=0:0.001:0.001*Nsteps;
figure;loglog(time(2:end),Ergodicity_Metric_save(2:end))
axis([0.001 5 0.0001,1])
xlabel('Time');
ylabel('Coverage Metric, \phi(t)');
title('Metric of ergodicity as a function of time')

