clear;close;clc
cd C:/Users/thmorin/Desktop/myModel/
addpath('./playground/');

playgroundConstants  %Loads up all the constants - should be functionalized
playgroundICs

T=293; %K
ft =1; %Temperature adjustment function
fm =1; %Moisture adjustment function
fph=1; %pH adjustment function. To be added later

aceV_1=nan(4,1);hV_1=nan(4,1); ch4V_1=nan(4,1);co2V_1=nan(4,1);acemu_1=nan(4,1);co2V_2=nan(4,1);h2V_2=nan(4,1);ch4V_2=nan(4,1);h2mu_2=nan(4,1);
fid=fopen('./out.dat','w');
% fprintf(fid,'Time,Acetate,CO2,CH4,H,M_CH4_ace\n');
fprintf(fid,'Time,Acetate,CO2,CH4,H2,H,M_CH4_ace,M_CH4_H2\n');
for t=1:c.nt-1
    for i=1:4
        if i==1; h=1; elseif i==2; h=0.5; elseif i==3; h=0.5; else; h=1;end
        %% Begin RK4 process - Generate fluxes here
%         [aceV_1(i),hV_1(i) ,ch4V_1(i),co2V_1(i),acemu_1(i)]=aceMethGens(m.CH4_ace.ace.Vmax,s.ace.tot(i),m.CH4_ace.Y,m.CH4_ace.ace.Km,ft,fm,fph);
        [aceV_1(i),hV_1(i),ch4V_1(i),co2V_1(i),acemu_1(i)]=aceMethGensTEST ...
            (m.CH4_ace.ace.Vmax,s.ace.tot(i),m.CH4_ace.ace.Km,ft,fm,fph,s.CH4.Gf,s.ace.Gf,...
                s.CO2.Gf,s.bio.Gf,c.m.Rm,c.m.Rw,s.CO2.tot(i),...
                s.CH4.tot(i),s.NH4.tot(i),s.PO3.tot(i),s.H.tot(i),c.R,T,c.CNrat,c.CPrat,m.CH4_ace.tot(i));

%         [co2V_2(i),h2V_2(i),ch4V_2(i),h2mu_2(i)           ]=h2MethGens (m.CH4_H2 .H2 .Vmax,s.H2.tot(i),m.CH4_H2.H2.Km,s.CO2.tot(i),m.CH4_H2.CO2.Km,m.CH4_H2.Y,ft,fm,fph);
%         fprintf('Writing shit %f and %f\n',co2V_2(i),h2V_2(i));
        %% Adjusting populations and concentrations
        s.ace.tot(i+1)=s.ace.tot(1) - c.dt*h*( aceV_1(i)*m.CH4_ace.tot(i) );
        s.CO2.tot(i+1)=s.CO2.tot(1) + c.dt*h*( co2V_1(i)*m.CH4_ace.tot(i) );%- co2V_2(i)*m.CH4_H2.tot(i) );
        s.CH4.tot(i+1)=s.CH4.tot(1) + c.dt*h*( ch4V_1(i)*m.CH4_ace.tot(i) );%+ ch4V_2(i)*m.CH4_H2.tot(i) );
%         s.H2 .tot(i+1)=s.H2 .tot(1) - c.dt*h*(  h2V_2(i)*m.CH4_H2. tot(i) );
        s.H  .tot(i+1)=s.H  .tot(1) - c.dt*h*(   hV_1(i)*m.CH4_ace.tot(1) );
        
        m.CH4_ace.tot(i+1) = m.CH4_ace.tot(1) + c.dt*h*m.CH4_ace.tot(i)*( acemu_1(i) - m.CH4_ace.Kd );
%         m.CH4_H2 .tot(i+1) = m.CH4_H2 .tot(1) + c.dt*h*m.CH4_H2 .tot(i)*(  h2mu_2(i) - m.CH4_H2 .Kd );
        s.NH4.tot(i+1)=s.NH4.tot(i);
        s.PO3.tot(i+1)=s.PO3.tot(i);
    end
    s.ace.tot(1)=1/6*( s.ace.tot(2) + 2*s.ace.tot(3) + 2*s.ace.tot(4) + s.ace.tot(5) );
    s.CO2.tot(1)=1/6*( s.CO2.tot(2) + 2*s.CO2.tot(3) + 2*s.CO2.tot(4) + s.CO2.tot(5) );
    s.CH4.tot(1)=1/6*( s.CH4.tot(2) + 2*s.CH4.tot(3) + 2*s.CH4.tot(4) + s.CH4.tot(5) );
%     s.H2 .tot(1)=1/6*( s.H2 .tot(2) + 2*s.H2 .tot(3) + 2*s.H2 .tot(4) + s.H2 .tot(5) );
    s.H  .tot(1)=1/6*( s.H  .tot(2) + 2*s.H  .tot(3) + 2*s.H  .tot(4) + s.H  .tot(5) );
    
    m.CH4_ace.tot(1)=1/6*( m.CH4_ace.tot(2) + 2*m.CH4_ace.tot(3) + 2*m.CH4_ace.tot(4) + m.CH4_ace.tot(5) );
%     m.CH4_H2 .tot(1)=1/6*( m.CH4_H2 .tot(2) + 2*m.CH4_H2 .tot(3) + 2*m.CH4_H2 .tot(4) + m.CH4_H2 .tot(5) );
    
%     fprintf(fid,'%f,%f,%f,%f,%f,%f\n',...
%         c.time(t),s.ace.tot(1),s.CO2.tot(1),s.CH4.tot(1),s.H.tot(1),m.CH4_ace.tot(1));
    fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f\n',...
        c.time(t),s.ace.tot(1),s.CO2.tot(1),s.CH4.tot(1),s.H2.tot(1),s.H.tot(1),m.CH4_ace.tot(1),m.CH4_H2.tot(1));
                      
%       fprintf('Time step %d written\n',c.time(t));
%     end
    if mod(t,100000)==0
        fprintf('Time step %d of %d written\n',c.time(t),c.nt-1);
    end
end
fclose(fid);

clear
data=importdata('out.dat');
for i=1:length(data.textdata)
    eval([data.textdata{:,i} '=data.data(:,i);']);
end

f1=figure();
f1.Units='Normalized';
f1.Position=[0 0 1 1];
subplot(2,1,1);
plot(Time,M_CH4_ace);
hold on;
plot(Time,M_CH4_H2);
grid on;xlabel('Time [d]');ylabel('Population');
legend('Acetoclastic methanogens','Hydrogenotrophic methanogens');

subplot(2,1,2);
plot(Time,Acetate); hold on;
plot(Time,CO2);
plot(Time,CH4);
plot(Time,H2);
grid on;xlabel('Time [d]');ylabel('Concentration');
legend('Acetate','CO_2','CH_4','H2');



        
    %     s.ace.tot(t+1)=s.ace.tot(t) - m.CH4_ace.tot(t)*m.CH4_ace.ace.V(t)*c.dt;    
    %     s.CO2.tot(t+1)=s.CO2.tot(t) + (m.CH4_ace.tot(t)*m.CH4_ace.CO2.V(t) - m.CH4_H2.tot(t)*m.CH4_H2.CO2.V(t))*c.dt;
    %     s.CH4.tot(t+1)=s.CH4.tot(t) + (m.CH4_ace.tot(t)*m.CH4_ace.CH4.V(t) + m.CH4_H2.tot(t)*m.CH4_H2.CH4.V(t))*c.dt;
%         s.H2 .tot(t+1)=s.H2 .tot(t) - m.CH4_H2.tot(t)*m.CH4_H2.H2.V(t)*c.dt;
%         s.H  .tot(t+1)=s.H  .tot(t) - m.CH4_ace.tot(t)*m.CH4_ace.H  .V(t)*c.dt;
%         m.CH4_ace.tot(t+1) = m.CH4_ace.tot(t)*(1 + (m.CH4_ace.mu(t) - m.CH4_ace.Kd)*c.dt );
%         m.CH4_H2 .tot(t+1) = m.CH4_H2 .tot(t)*(1 + (m.CH4_H2 .mu(t) - m.CH4_H2 .Kd)*c.dt );