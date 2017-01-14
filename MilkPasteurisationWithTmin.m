clear variables; close all; clc;

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

delta_PAST5=linspace(0.01,10,9999);
T_crpast_c=T_cream_35+delta_PAST5;
T_crpast_a=T_crpast_b-delta_PAST5;
m_crpast_c=m_cream_35;
m_crpast_a=m_cream_35;
cp_crpast_c=cp_cream_35;
cp_crpast_a=cp_cream_35;

%Evaluation of extra cost for heat exchanger modification

k1_ex=3.8528;
k2_ex=0.4242;
Q_PAST5=m_cream_35*cp_cream_35*(T_crpast_a-T_cream_35)*1000;

LMTD_PAST5=(T_crpast_b-T_crpast_a);
alpha=1000;
U=alpha/2;
A_PAST5=Q_PAST5./(U*LMTD_PAST5);

It=556.8;
Iref=389.5;
FBM=4.74;
e=0.96;
i=0.08;
n=20;

Cp_PAST5=It/Iref*10.^(k1_ex+k2_ex*log10(A_PAST5));
CBM_PAST5=FBM*Cp_PAST5*e;
IC_PAST5=CBM_PAST5*i*(1+i)^n/((1+i)^n-1);


%-

T_cream=4;
m_cream=m_cream_35;
cp_cream=cp_cream_35;

%Supposition 2:

T_past_b=75;
m_past_b=m_milk_0;
cp_past_b=cp_milk_0;

delta_PAST2=linspace(0.01,10,9999);
T_past_c=T_milk_0+delta_PAST2;
T_past_a=T_past_b-delta_PAST2;
m_past_c=m_milk_0;
m_past_a=m_milk_0;
cp_past_c=cp_milk_0;
cp_past_a=cp_milk_0;
T_past_d=T_past_c-(m_ref_in*cp_ref_in/cp_past_c/m_past_c)*(T_past_cen-T_ref_out);

%Evaluation of extra cost for heat exchanger modification

k1_ex=3.8528;
k2_ex=0.4242;
Q_PAST2=m_milk_0*cp_milk_0*(-T_milk_0+T_past_a)*1000;

LMTD_PAST2=(T_past_b-T_past_a);
alpha=1000;
U=alpha/2;
A_PAST2=Q_PAST2./(U*LMTD_PAST2);

Q_PAST1=m_ref_out*cp_ref_in*(-T_ref_out+T_past_cen)*1000;
LMTD_PAST1=(T_past_c-T_past_cen-T_past_d+T_ref_out)./log((T_past_c-T_past_cen)./(T_past_d-T_ref_out));
A_PAST1=Q_PAST1./(U*LMTD_PAST1);

It=556.8;
Iref=389.5;
FBM=4.74;
e=0.96;
i=0.08;
n=20;

Cp_PAST2=It/Iref*10.^(k1_ex+k2_ex*log10(A_PAST2));
Cp_PAST1=It/Iref*10.^(k1_ex+k2_ex*log10(A_PAST1));
CBM_PAST2=FBM*Cp_PAST2*e;
CBM_PAST1=FBM*Cp_PAST1*e;
IC_PAST2=CBM_PAST2*i*(1+i)^n/((1+i)^n-1);
IC_PAST1=CBM_PAST1*i*(1+i)^n/((1+i)^n-1);

delta_PAST1a=(T_past_c-T_past_cen);
delta_PAST1b=(T_past_d-T_ref_out);
%-

%T_past_d=T_past_c-(m_ref_in*cp_ref_in/cp_past_c/m_past_c)*(T_past_cen-T_ref_out);
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

valcp = cp_crpast_a
valm = m_crpast_a
valTb = T_crpast_b
valTa = T_crpast_a(1:10)
valOT = ot
valgp = gp



%Optimal deltaT

OCbase=REF;
OC1=PAST6 + PAST7;
OC2=PAST4 + PAST3 + OCbase;

[minOC1, index1]=min(OC1+IC_PAST5);
[minOC2, index2]=min(OC2+IC_PAST2+IC_PAST1);

deltaT_PAST2min=delta_PAST2(index2);
deltaT_PAST5min=delta_PAST5(index1);
deltaT_PAST1min=min(delta_PAST1a(index2),delta_PAST1b(index2));

Total=OCbase + PAST6(index1) + PAST7(index1) + PAST4(index2) + PAST3(index2);


%%
t_payback_1 = CBM_PAST5(index1)/(OC1(9999)-OC1(index1));
t_payback_2 = (CBM_PAST2(index2)+CBM_PAST1(index2))/(OC2(9999)-OC2(index2));

t_payback_1_2 = (CBM_PAST1(index2)+CBM_PAST2(index2)+CBM_PAST5(index1))/(OC1(9999)+OC2(9999)-OC1(index1)-OC2(index2));

figure(1);
hold on;
plot(delta_PAST5,(IC_PAST5+OC1)/1000,'g');
plot(delta_PAST5, IC_PAST5/1000,'r');
plot(delta_PAST5, OC1/1000,'b');
plot([deltaT_PAST5min deltaT_PAST5min],[0 max((IC_PAST5+OC1)/1000)],'--k')
legend('Total cost','Investment cost','Operating cost');
xlabel('\Delta T_{min} [°C]');
ylabel('Cost per year [kCHF/year]');
title('Part 1');    

figure(2);
hold on;
plot(delta_PAST2,(IC_PAST2+IC_PAST1+OC2)/1000,'g');
plot(delta_PAST2, (IC_PAST2+IC_PAST1)/1000,'r');
plot(delta_PAST2, OC2/1000,'b');
plot([deltaT_PAST2min deltaT_PAST2min],[0 max((IC_PAST2+IC_PAST1+OC2)/1000)],'--k')
legend('Total cost','Investment cost','Operating cost');
xlabel('\Delta T_{min} [°C]');
ylabel('Cost per year [kCHF/year]');
title('Part 2');

figure(3);
hold on;
x=3192;

plot(delta_PAST2(x:1:9999),(IC_PAST2(x:1:9999)+IC_PAST1(x:1:9999)+OC2(x:1:9999))/1000,'g');
plot(delta_PAST2(x:1:9999), (IC_PAST2(x:1:9999)+IC_PAST1(x:1:9999))/1000,'r');
plot(delta_PAST2(x:1:9999), OC2(x:1:9999)/1000,'b');
plot([deltaT_PAST2min deltaT_PAST2min],[0 max((IC_PAST2+IC_PAST1+OC2)/1000)],'--k')
legend('Total cost','Investment cost','Operating cost');
xlabel('\Delta T_{min} [°C]');
ylabel('Cost per year [kCHF/year]');
title('Part 2 limite');

