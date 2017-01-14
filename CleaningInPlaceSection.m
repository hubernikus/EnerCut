clear all;
close all;
clc;

% 1) Mass&Energy balance :

T_water=25;

T_cip_e1=70;
m_cip_e1=4;
cp_water_70 = 4.191;
cp_water_25 = 4.180;
cp_water_75 = 4.194;
cp_water_80 = 4.198;

% If S-5 is isothermal, and cp==constant then :

T_recup_2=75;

syms m_water_var m_recup_2_var

[m_water m_recup_2]=solve(m_water_var*T_water + m_recup_2_var*T_recup_2 == m_cip_e1*T_cip_e1 , m_water_var + m_recup_2_var == m_cip_e1 , m_water_var, m_recup_2_var);

% 2) Energy bill : 

ot=8000;
ep=0.081;
gp=0.038;

Steam=(m_cip_e1*cp_water_75*(80-70))*ot*gp/0.9;

Total=Steam;
