%-------------------------------------------------------------------------------
%---------------------------PARAMETERS------------------------------------------
%FIXME: Switch this to a reader from an input sheet at some point
%                       Simulation domain
c.dt=1;
c.nCol=1;     c.nRow=1;     c.nLyr=1; %For when I expand it later, form should be m(row,col,lyr).XXX.YYY etc
c.R=8.29; %J mol^-1 K^-1

%                     Equilibration constants
p.NO3.kH=1;   p.NO2.kH=1;   p.NH4.kH=1;   p.CH4.kH=1;   p.CO2.kH=1;   p.but.kH=1;
p.mon.kH=1;   p.ace.kH=1;   p.H2.kH=1;    p.O2.kH=1;    p.PO3.kH=1;   p.PO2.kH=1;
p.PO4.kH=1;   p.eth.kH=1;   p.meth.kH=1;
p.micC.kH=1;  p.micN.kH=1;  p.micP.kH=1;
p.names=fieldnames(p);

%                       Molecular weights (g mol^-1)
p.NO3.MW=62;  p.NO2.MW=46;  p.NH4.MW=18;  p.CH4.MW=16;  p.CO2.MW=44;
p.O2.MW=32;   p.ace.MW=;    p.H2.MW =2;   p.O2.MW =32;  p.PO3.MW=;  p.PO2.MW=;
p.PO4.MW=;  p.eth.MW=;  p.meth.MW=;
p.micC.MW=;  p.micN.MW;  p.micP.MW;
%-------------------------------------------------------------------------------