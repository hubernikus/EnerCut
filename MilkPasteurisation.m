clear all;
close all;
clc;

% 1) Mass&Energy balance :

T_ref_in=6; 
m_ref_in=8;

T_ref_out=4;
m_ref_out=m_ref_in;

T_past_cen=60;
m_past_cen=m_ref_out;

T_cream_35=T_past_cen;
m_cream_35=0.48;
cp_cream_35=3.4;

T_milk_0=T_past_cen;
m_milk_0=7.52;
cp_milk_0=3.8;

cp_ref_in=(m_cream_35*cp_cream_35+m_milk_0*cp_milk_0)/m_ref_in;

%Supposition 1:

T_crpast_b=75;
m_crpast_b=m_cream_35;
cp_crpast_b=cp_cream_35;

delta_PAST5=10;
T_crpast_c=T_cream_35+delta_PAST5;
T_crpast_a=T_crpast_b-delta_PAST5;
m_crpast_c=m_cream_35;
m_crpast_a=m_cream_35;
cp_crpast_c=cp_cream_35;
cp_crpast_a=cp_cream_35;

%-

T_cream=4;
m_cream=m_cream_35;
cp_cream=cp_cream_35;

%Supposition 2:

T_past_b=75;
m_past_b=m_milk_0;
cp_past_b=cp_milk_0;

delta_PAST2=10;
T_past_c=T_milk_0+delta_PAST2;
T_past_a=T_past_b-delta_PAST2;
m_past_c=m_milk_0;
m_past_a=m_milk_0;
cp_past_c=cp_milk_0;
cp_past_a=cp_milk_0;

%-

T_past_d=T_past_c-(m_ref_in*cp_ref_in/cp_past_c/m_past_c)*(T_past_cen-T_ref_out);
m_past_d=m_milk_0;
cp_past_d=cp_milk_0;

T_milk=4;
m_milk=m_milk_0;
cp_milk=cp_milk_0;

% 2) Current energy bill

COP=4;
ot=8000;
ep=0.081;
gp=0.038;

REF=m_ref_in*cp_ref_in*(T_ref_in-T_ref_out)/COP*ot*ep;
PAST7=cp_cream*m_cream*(T_crpast_c-T_cream)/COP*ot*ep;
PAST4=cp_milk*m_milk*(T_past_d-T_milk)/COP*ot*ep;

PAST3=cp_past_a*m_past_a*(T_past_b-T_past_a)*ot*gp/0.9;
PAST6=cp_crpast_a*m_crpast_a*(T_crpast_b-T_crpast_a)*ot*gp/0.9;

Total=REF + PAST7 + PAST4 + PAST3 + PAST6;

C1=PAST7 + PAST6;
C2=PAST3 + PAST4 + REF;





