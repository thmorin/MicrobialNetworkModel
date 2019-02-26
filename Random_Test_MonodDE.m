%   Extensive testing of monod diff eq solver by generating random values
%   of parameters over an estimated range of valeus 

%Defining range of test Units: days
range = [0 200];

%   initial conditions
s  = 500;   %   Substrate concentration
m1 = 10;    %   Microbial population 1
n  = 10000;   %   Nitrogen concentration
p  = 100;   %   Phosphorus concentration
a  = 100;     %   Acetate concentration
m2 = 15;    %   Microbial population 2
m3 = 5;
CH4= 0;
%   Combine into single vector for input into solver
IC = [s;m1;n;p;a;m2;m3;CH4];

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

%   For loop to generate random parameters for each microbe population goes
%   from 1 to 2 b/c only set up for 2 populations at the moment
for i = 1:3
    %   Randomly generating parameters
   ks(i) = ksmin+rand(1,1)*(ksmax-ksmin); 
   kd(i) = kdmin+rand(1,1)*(kdmax-kdmin);
   mu(i) = mumin+rand(1,1)*(mumax-mumin);
   yield(i) = ymin+rand(1,1)*(ymax-ymin);
   kn(i) = knmin+rand(1,1)*(knmax-knmin);
   kp(i) = kpmin+rand(1,1)*(kpmax-kpmin);
    %   Populating 'params' with randomly generated parameters
   params(i) = parameters(ks(i),kd(i),mu(i),yield(i),kn(i),kp(i))
end 

%   Testing different ode solvers
%[tv,Yv]=ode15s(@(tv,Yv) monodDE(tv,Yv,params),range,IC);
%[tv,Yv]=ode23(@(tv,Yv) monodDE(tv,Yv,params),range,IC);
[tv,Yv]=ode23s(@(tv,Yv) monodDE(tv,Yv,params),range,IC);
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
plot(tv, Yv(:,2))
hold on
plot(tv,Yv(:,6),'-r')
plot(tv,Yv(:,7),'-c');
title('Microbial Population vs Time')
xlabel('Time (days)')
ylabel('Microbial Population')

subplot(6,1,3)
plot(tv,Yv(:,3))
title('Nitrogen Concentration vs Time')
xlabel('Time (days)')
ylabel('Nitrogen Concentration')

subplot(6,1,4)
plot(tv,Yv(:,4))
title('Phosphorus Concentration vs Time')
xlabel('Time (days)')
ylabel('Phosphorus Concentration')

subplot(6,1,5)
plot(tv,Yv(:,5))
title('Acetate Concentration vs Time')
xlabel('Time (days)')
ylabel('Acetate Conc.')

subplot(6,1,6)
plot(tv,Yv(:,8))
title('CH_4 Concentration vs Time')
xlabel('Time (days)')
ylabel('CH_4 Conc.')

params(1)
params(2)