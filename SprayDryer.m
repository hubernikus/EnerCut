clear all;
close all;
clc;

% 1) Energy&Mass balance :

%Air at ambient temperature of 20°C

T_air=20;
cp_air_20=1.005;
cp_air_100=1.009;
p_air=10^5;

T_concentrated_milk=4;
cp_concentrated_milk=2.61;
m_concentrated_milk=1.21;

T_milk_powder=70;
cp_milk_powder=1.165;
m_milk_powder=m_concentrated_milk*0.5/0.96;

cp_moyen_concentrated_powder=(cp_milk_powder+cp_concentrated_milk)/2;

m_water=m_concentrated_milk-m_milk_powder;
cp_water_20=4.182;

m_air_sec=(-cp_milk_powder*m_milk_powder*(70+273)-cp_water_20*m_water*(20+273)+cp_concentrated_milk*m_concentrated_milk*(4+273))/(cp_air_100*(200-20));

m_air=m_air_sec/0.98;


% 2) Energy bill

COP=4;
ot=8000;
ep=0.081;
gp=0.038;

Glycol_cooling=(m_air*cp_air_20*(20-4))/COP*ot*ep;

Boiler=(m_air_sec*cp_air_100*(200-4))*ot*gp/0.9;

Total=Boiler + Glycol_cooling;



