clear;close;clc
%   Extensive testing of monod diff eq solver by generating random values
%   of parameters over an estimated range of valeus 

%Defining range of test Units: days
range = [0 200];

%   initial conditions
S  = [500;0;0;0;0];%500*ones(5,1);   %   Substrate concentration
N  = 100000*ones(5,1);   %   Nitrogen concentration
P  = 1*ones(5,1);   %   Phosphorus concentration
A  = 100*ones(5,1);     %   Acetate concentration
CH4= 0*ones(5,1);

% Replicate Zerkle Profiles coming up with new model grid
S2O3=[0;0.1;]


M1 = 10*ones(5,1);    %   Microbial population 1
M2 = 15*ones(5,1);    %   Microbial population 2
M3 = 5*ones(5,1);

%   Combine into single vector for input into solver
IC = [S N P A CH4 M1 M2 M3];

%   Defining Ks Range
ksmin = 1;
ksmax  = 500;
%   Defining Kd Range
kdmin = .0001;
kdmax = .25;
%   Defining max growth rate range
mumax = 25;
mumin = 1;
%   Defining Yield factor
ymin = .0001;
ymax = .5;
%   Defining kN
knmin = 1;
knmax = 200;
%   Defining kP
kpmin = 1;
kpmax = 100;

%   Defining Diffusion parameter
Dmin = 1;
Dmax = 100;

%   For loop to generate random parameters for each microbe population goes
%   from 1 to 2 b/c only set up for 2 populations at the moment
for i = 1:3
    %   Randomly generating parameters
   ks(i)    = ksmin + rand(1,1)*(ksmax-ksmin); 
   kd(i)    = kdmin + rand(1,1)*(kdmax-kdmin);
   mu(i)    = mumin + rand(1,1)*(mumax-mumin);
   yield(i) = ymin  + rand(1,1)*(ymax-ymin);
   kn(i)    = knmin + rand(1,1)*(knmax-knmin);
   kp(i)    = kpmin + rand(1,1)*(kpmax-kpmin);
    %   Populating 'params' with randomly generated parameters
   m(i) = parameters(ks(i),kd(i),mu(i),yield(i),kn(i),kp(i));
end 

for i = 1:5
    s(i).D = Dmin+rand(1,1)*(Dmax-Dmin);  %1=DOC
    s(i).Gf = -1026.55 %Phosphate kJ/mol
    
end

initialConds
p.Z=[5;10;15;20;25;30];
p.dz=diff(p.Z);

SO2_IC=interp1(data(:,1),data(:,2),p.Z);



%   Testing different ode solvers
%[tv,Yv]=ode15s(@(tv,Yv) monodDE(tv,Yv,params),range,IC);
%[tv,Yv]=ode23(@(tv,Yv) monodDE(tv,Yv,params),range,IC);
[tv,Yv]=ode23s(@(tv,Yv) Tim_monodDE(tv,Yv,m,s,p),range,IC);
%[tv,Yv]=ode45(@(tv,Yv) monodDE(tv,Yv,params),range,IC);

%   Plotting results from ode solver
f1=figure();
f1.Units='Normalized';
f1.Position=[0 0 1 1];
subplot(6,1,1)
plot(tv, Yv(:,1))
title('Substrate Concentration vs Time')
xlabel('Time (days)')
ylabel('Substrate Concentration')

subplot(6,1,2)
plot(tv, Yv(:,(6-1)*5+1))
hold on
plot(tv,Yv(:,(7-1)*5+1),'-r')
plot(tv,Yv(:,(8-1)*5+1),'-c');
title('Microbial Population vs Time')
xlabel('Time (days)')
ylabel('Microbial Population')

subplot(6,1,3)
plot(tv,Yv(:,(2-1)*5+1))
title('Nitrogen Concentration vs Time')
xlabel('Time (days)')
ylabel('Nitrogen Concentration')

subplot(6,1,4)
plot(tv,Yv(:,(3-1)*5+1))
title('Phosphorus Concentration vs Time')
xlabel('Time (days)')
ylabel('Phosphorus Concentration')

subplot(6,1,5)
plot(tv,Yv(:,(4-1)*5+1))
title('Acetate Concentration vs Time')
xlabel('Time (days)')
ylabel('Acetate Conc.')

subplot(6,1,6)
plot(tv,Yv(:,(5-1)*5+1))
title('CH_4 Concentration vs Time')
xlabel('Time (days)')
ylabel('CH_4 Conc.')

figure();
image(tv,1:5,Yv(:,1:5)');colorbar();
xlabel('Time (days)');
ylabel('Depth');