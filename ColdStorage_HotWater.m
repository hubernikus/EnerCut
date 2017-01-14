clear all;
close all;
clc;

% 1) Energy&Mass balance : 

% a) Cold Storage :

Q = 500*10^3;

% b) Hot Water : 

cp_hw_moy = 4.178;
T_hw_in=15;
T_hw_out=55;
m_hw_in=1;
m_hw_out=1;

% 2) Energy bill : 

COP=4;
ot=8000;
ep=0.081;
gp=0.038;

ColdStorage = 500/COP*ot*ep;

HotWater = m_hw_in*cp_hw_moy*(T_hw_out-T_hw_in)*ot*gp/0.9;

Total = HotWater + ColdStorage; 