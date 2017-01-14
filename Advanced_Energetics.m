clear all;
close all;
clc;

disp 'Starting  <<Advanced Energetics - Project>> ';

%% Pasteurisation section

disp 'Starting  <<Pasteurisation Section>> ';
for Country=1:2;

for ProcessType=1:2;

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


Q_PAST2=m_milk_0*cp_milk_0*(-T_milk_0+T_past_a)*1000;

LMTD_PAST2=(T_past_b-T_past_a);
alpha=1000;
U=alpha/2;
A_PAST2=Q_PAST2./(U*LMTD_PAST2);

Q_PAST1=m_ref_out*cp_ref_in*(-T_ref_out+T_past_cen)*1000;
LMTD_PAST1=(T_past_c-T_past_cen-T_past_d+T_ref_out)./log((T_past_c-T_past_cen)./(T_past_d-T_ref_out));
A_PAST1=Q_PAST1./(U*LMTD_PAST1);

Cp_PAST2=It/Iref*10.^(k1_ex+k2_ex*log10(A_PAST2));
Cp_PAST1=It/Iref*10.^(k1_ex+k2_ex*log10(A_PAST1));
CBM_PAST2=FBM*Cp_PAST2*e;
CBM_PAST1=FBM*Cp_PAST1*e;
IC_PAST2=CBM_PAST2*i*(1+i)^n/((1+i)^n-1);
IC_PAST1=CBM_PAST1*i*(1+i)^n/((1+i)^n-1);

delta_PAST1a=(T_past_c-T_past_cen);
delta_PAST1b=(T_past_d-T_ref_out);
%-

m_past_d=m_milk_0;
cp_past_d=cp_milk_0;

T_milk=4;
m_milk=m_milk_0;
cp_milk=cp_milk_0;

% 2) Current energy bill

COP=4;
ot=8000;
ep=[0.081 0.124];
gp=[0.038 0.05];

if Country==1;
    ep=ep(1);
    gp=gp(1);
end
if Country==2;
    ep=ep(2);
    gp=gp(2);
end

REF=m_ref_in*cp_ref_in*(T_ref_in-T_ref_out)/COP*ot*ep;
PAST7=cp_cream*m_cream*(T_crpast_c-T_cream)/COP*ot*ep;
PAST4=cp_milk*m_milk*(T_past_d-T_milk)/COP*ot*ep;

PAST3=cp_past_a*m_past_a*(T_past_b-T_past_a)*ot*gp/0.9;
PAST6=cp_crpast_a*m_crpast_a*(T_crpast_b-T_crpast_a)*ot*gp/0.9;



%Optimal deltaT

OCbase=REF;
OC1=PAST6 + PAST7;
OC2=PAST4 + PAST3 + OCbase;

[minOC1, index1]=min(OC1+IC_PAST5);
[minOC2, index2]=min(OC2+IC_PAST2+IC_PAST1);

if ProcessType==1;
    index1=9999;
    index2=9999;
end 

deltaT_PAST2min=delta_PAST2(index2);
deltaT_PAST5min=delta_PAST5(index1);
deltaT_PAST1min=min(delta_PAST1a(index2),delta_PAST1b(index2));

Total=OCbase + PAST6(index1) + PAST7(index1) + PAST4(index2) + PAST3(index2);

t_payback_1 = CBM_PAST5(index1)/(OC1(9999)-OC1(index1));
t_payback_2 = (CBM_PAST2(index2)+CBM_PAST1(index2))/(OC2(9999)-OC2(index2));

t_payback_1_2 = (CBM_PAST1(index2)+CBM_PAST2(index2)+CBM_PAST5(index1))/(OC1(9999)+OC2(9999)-OC1(index1)-OC2(index2));


%Energy balances:

Q_REF=cp_ref_in*m_ref_in*(T_ref_in-T_ref_out);
Q_PAST1=cp_milk*m_milk*(T_past_c-T_past_d);
Q_PAST1=Q_PAST1(index2);
Q_PAST2=cp_milk*m_milk*(T_past_a-T_milk_0);
Q_PAST2=Q_PAST2(index2);
Q_PAST3=cp_milk*m_milk*(T_past_b-T_past_a); 
Q_PAST3=Q_PAST3(index2);
Q_PAST4=cp_milk*m_milk*(T_past_d-T_milk);
Q_PAST4=Q_PAST4(index2);
Q_PAST5=cp_cream_35*m_cream_35*(T_crpast_a-T_cream_35);
Q_PAST5=Q_PAST5(index1);
Q_PAST6=cp_cream_35*m_cream_35*(T_crpast_b-T_crpast_a);
Q_PAST6=Q_PAST6(index1);
Q_PAST7=cp_cream_35*m_cream_35*(T_crpast_c-T_cream);
Q_PAST7=Q_PAST7(index1);

close all;

if ProcessType==2;
figure(1);
hold on;
plot(delta_PAST5,(IC_PAST5+OC1)/1000,'g');
plot(delta_PAST5, IC_PAST5/1000,'r');
plot(delta_PAST5, OC1/1000,'b');
plot([deltaT_PAST5min deltaT_PAST5min],[0 max((IC_PAST5+OC1)/1000)],'--k');
legend('Total cost','Investment cost','Operating cost');
xlabel('\Delta T_{min} [°C]');
ylabel('Cost per year [kCHF/year]');
if Country==1;
title('Cream Process Optimization France');
print('Cream_Process_Optimization_France','-dpng');
end
if Country==2;
title('Cream Process Optimization Germany');
print('Cream_Process_Optimization_Germany','-dpng');
end

figure(2);
hold on;
plot(delta_PAST2,(IC_PAST2+IC_PAST1+OC2)/1000,'g');
plot(delta_PAST2, (IC_PAST2+IC_PAST1)/1000,'r');
plot(delta_PAST2, OC2/1000,'b');
plot([deltaT_PAST2min deltaT_PAST2min],[0 max((IC_PAST2+IC_PAST1+OC2)/1000)],'--k');
legend('Total cost','Investment cost','Operating cost');
xlabel('\Delta T_{min} [°C]');
ylabel('Cost per year [kCHF/year]');
if Country==1;
title('Milk Process Optimization France');
print('Milk_Process_Optimization_France','-dpng');
end
if Country==2;
title('Milk Process Optimization Germany');
print('Milk_Process_Optimization_Germany','-dpng');
end

figure(3);
hold on;
%x=3192;


for RechercheMinimum=1:9999;
    if T_past_d(RechercheMinimum)-T_ref_out>0;
        x=RechercheMinimum;
        break
    end
end


plot(delta_PAST2(x:1:9999),(IC_PAST2(x:1:9999)+IC_PAST1(x:1:9999)+OC2(x:1:9999))/1000,'g');
plot(delta_PAST2(x:1:9999), (IC_PAST2(x:1:9999)+IC_PAST1(x:1:9999))/1000,'r');
plot(delta_PAST2(x:1:9999), OC2(x:1:9999)/1000,'b');
plot([deltaT_PAST2min deltaT_PAST2min],[0 max((IC_PAST2+IC_PAST1+OC2)/1000)],'--k');
legend('Total cost','Investment cost','Operating cost');
xlabel('\Delta T_{min} [°C]');
ylabel('Cost per year [kCHF/year]');
if Country==1;
title('Milk Process Optimization Corrected France');
print('Milk_Process_Optimization_Corrected_France','-dpng');
end
if Country==2;
title('Milk Process Optimization Corrected Germany');
print('Milk_Process_Optimization_Corrected_Germany','-dpng');
end

close all;

end


S1.Name = {'ref_in'; 'ref_out'; 'past_cen'; 'milk_0'; 'past_a'; 'past_b'; 'past_c'; 'past_d'; 'milk'; 'cream_35'; 'crpast_a'; 'crpast_b'; 'crpast_c'; 'cream'};
S1.MassFlow = [m_ref_in; m_ref_in; m_ref_in;m_milk_0;m_milk_0;m_milk_0;m_milk_0;m_milk_0;m_milk_0;m_cream_35;m_cream_35;m_cream_35;m_cream_35;m_cream_35;];
S1.Temperature = [T_ref_in;T_ref_out;T_past_cen;T_milk_0;T_past_a(index2);T_past_b;T_past_c(index2);T_past_d(index2);T_milk;T_cream_35;T_crpast_a(index1);T_crpast_b;T_crpast_c(index1);T_cream];
S1.Cp = [cp_ref_in;cp_ref_in;cp_ref_in;cp_milk_0;cp_milk_0;cp_milk_0;cp_milk_0;cp_milk_0;cp_milk_0;cp_cream_35;cp_cream_35;cp_cream_35;cp_cream_35;cp_cream_35;];



if ProcessType==1;
Mass_Temperature_10 = struct2table(S1);
writetable(Mass_Temperature_10,'Mass_Temperature_10.xls');

S2.Name = {'REF';'PAST1';'PAST2';'PAST3';'PAST4';'PAST5';'PAST6';'PAST7'};
S2.HeatTransfer = [Q_REF;Q_PAST1;Q_PAST2;Q_PAST3;Q_PAST4;Q_PAST5;Q_PAST6;Q_PAST7];

Heat_Transfer_10 = struct2table(S2);
writetable(Heat_Transfer_10,'Heat_Transfer_10.xls');
end


if ProcessType==2;
    
    if Country==1;
Mass_Temperature = struct2table(S1);
writetable(Mass_Temperature,'Mass_Temperature_France.xls');

S2.Name = {'REF';'PAST1';'PAST2';'PAST3';'PAST4';'PAST5';'PAST6';'PAST7'};
S2.HeatTransfer = [Q_REF;Q_PAST1;Q_PAST2;Q_PAST3;Q_PAST4;Q_PAST5;Q_PAST6;Q_PAST7];

Heat_Transfer = struct2table(S2);
writetable(Heat_Transfer,'Heat_Transfer_France.xls');
    end
    if Country==2;
Mass_Temperature = struct2table(S1);
writetable(Mass_Temperature,'Mass_Temperature_Germany.xls');

S2.Name = {'REF';'PAST1';'PAST2';'PAST3';'PAST4';'PAST5';'PAST6';'PAST7'};
S2.HeatTransfer = [Q_REF;Q_PAST1;Q_PAST2;Q_PAST3;Q_PAST4;Q_PAST5;Q_PAST6;Q_PAST7];

Heat_Transfer = struct2table(S2);
writetable(Heat_Transfer,'Heat_Transfer_Germany.xls');
    end
end




if ProcessType==2;

    S3.Unit={'REF';'PAST1';'PAST2';'PAST3';'PAST4';'PAST5';'PAST6';'PAST7';'Total'};
    S3.Cost=[REF;0;0;PAST3(index2);PAST4(index2);0;PAST6(index1);PAST7(index1);REF + PAST6(index1) + PAST7(index1) + PAST4(index2) + PAST3(index2)];
    S3.Payback={'/';'/';'/';'/';'/';'/';'/';'/';num2str(t_payback_1_2)};
    S3.DeltaTmin={'/';num2str(deltaT_PAST1min);num2str(deltaT_PAST2min);'/';'/';num2str(deltaT_PAST5min);'/';'/';'/'};
    Economics = struct2table(S3);
    
if Country==1;
   
     Cost_Pasteurisation_Fr_Optimized=REF + PAST6(index1) + PAST7(index1) + PAST4(index2) + PAST3(index2);
    writetable(Economics,'Cost_Analysis_France_Optimized.xls'); 
        
end
  if Country==2;
    writetable(Economics,'Cost_Analysis_Germany_Optimized.xls');  
    Cost_Pasteurisation_Ge_Optimized=REF + PAST6(index1) + PAST7(index1) + PAST4(index2) + PAST3(index2);
end

end

 if ProcessType==1;
     
      S3.Unit={'REF';'PAST1';'PAST2';'PAST3';'PAST4';'PAST5';'PAST6';'PAST7';'Total'};
    S3.Cost=[REF;0;0;PAST3(9999);PAST4(9999);0;PAST6(9999);PAST7(9999);REF + PAST6(9999) + PAST7(9999) + PAST4(9999) + PAST3(9999)];
    S3.DeltaTmin={'/';num2str(deltaT_PAST1min);num2str(deltaT_PAST2min);'/';'/';num2str(deltaT_PAST5min);'/';'/';'/'};
    Economics = struct2table(S3);
    
    if Country==1;
   
    writetable(Economics,'Cost_Analysis_France_Normal.xls'); 
     Cost_Pasteurisation_Fr= REF + PAST6(9999) + PAST7(9999) + PAST4(9999) + PAST3(9999) ; 
end
  if Country==2;
    writetable(Economics,'Cost_Analysis_Germany_Normal.xls');  
     Cost_Pasteurisation_Ge= REF + PAST6(9999) + PAST7(9999) + PAST4(9999) + PAST3(9999) ; 
end
 end

    end
   close all; 
end

disp 'Ending  <<Pasteurisation section>> ';
%% Spray Dryer Section

% 1) Energy&Mass balance :
disp 'Starting  <<Dryer section>> ';
for DryerCountry=(1:2);
    

%Air at ambient temperature of 20°C

T_air=25;
cp_air_20=1.005;
cp_air_100=1.009;
p_air=10^5;
cp_air=1.005;

T_concentrated_milk=4;
cp_concentrated_milk=2.61;
m_concentrated_milk=1.21;

T_milk_powder=70;
cp_milk_powder=1.165;
m_milk_powder=m_concentrated_milk*0.5/0.96;

cp_moyen_concentrated_powder=(cp_milk_powder+cp_concentrated_milk)/2;

m_water=m_concentrated_milk-m_milk_powder;
cp_water_20=4.182;
cp_water=4.182;

T_a_exit=25;
T_water_exit=25;
%m_air_sec=(cp_milk_powder*m_milk_powder*(70+273)+cp_water_20*m_water*(20+273)-cp_concentrated_milk*m_concentrated_milk*(4+273)+cp_air_20*(20+273))/(cp_air_100*(200+273));

%m_air=m_air_sec/0.98;

T_air_exit=(4:0.01:87);

Pv=exp(-3928.5./(231.667+T_air_exit))*140975;

w=(18/28.9)*Pv./(1.013-Pv);

m_air_sec=m_water./w;

E=cp_air*m_air_sec*(200+273)+cp_concentrated_milk*4*m_concentrated_milk-cp_air*m_air_sec.*(T_air_exit+273)-T_air_exit*cp_water*m_water-2257*m_water-cp_milk_powder*m_milk_powder*(70+273);

for Recherche0=1:length(T_air_exit);
    if E(Recherche0)<0;
        x2=Recherche0-1;
        break
    end
end

m_air_sec=m_air_sec(x2);
T_air_exit=T_air_exit(x2);
Pv=Pv(x2);
w=w(x2);

m_air=m_air_sec/0.98;

m_water_separator=0.02*m_air;
T_separator=4;

% 2) Energy bill

COP=4;
ot=8000;
ep=[0.081 0.124];
gp=[0.038 0.05];

Q_GycolCooling=cp_air_20*m_air*(25-4)+2257*m_air*0.02;
Q_Heating_by_NG=cp_air_100*m_air_sec*(200-4);

if DryerCountry==1;

Glycol_cooling=(m_air*cp_air_20*(25-4)+2257*m_air*0.02)/COP*ot*ep(1);

Boiler=(m_air_sec*cp_air_100*(200-4))*ot*gp(1)/0.9;

Total=Boiler + Glycol_cooling;

    S4.Unit={'GlycolCooling';'HeatingByNG';'Total';'PumpedAir';'ConcentratedMilk';'MilkPowder';'DryAir';'ExtractedWater'};
    S4.Cost={num2str(Glycol_cooling);num2str(Boiler);num2str(Total);'/';'/';'/';'/';'/'};
    S4.HeatTransfer={num2str(Q_GycolCooling);num2str(Q_Heating_by_NG);'/';'/';'/';'/';'/';'/'};
    S4.Temperature={'/';'/';'/';num2str(T_air);num2str(T_concentrated_milk);num2str(T_milk_powder);num2str(T_air_exit);num2str(T_air_exit)};
    S4.MassFlow={'/';'/';'/';num2str(m_air);num2str(m_concentrated_milk);num2str(m_milk_powder);num2str(m_air_sec);num2str(m_water)};
    S4.Cp={'/';'/';'/';num2str(cp_air_20);num2str(cp_concentrated_milk);num2str(cp_milk_powder);num2str(cp_air_20);num2str(cp_water_20)};
    
    DryerValues = struct2table(S4);

    writetable( DryerValues,'Spray_Dryer_Values_France.xls');
    Cost_Spray_dryer_Fr=double(Total);
end

if DryerCountry==2;

Glycol_cooling=(m_air*cp_air_20*(25-4))/COP*ot*ep(2);

Boiler=(m_air_sec*cp_air_100*(200-4))*ot*gp(2)/0.9;

Total=Boiler + Glycol_cooling;

    S4.Unit={'GlycolCooling';'HeatingByNG';'Total';'PumpedAir';'ConcentratedMilk';'MilkPowder';'DryAir';'ExtractedWater'};
    S4.Cost={num2str(Glycol_cooling);num2str(Boiler);num2str(Total);'/';'/';'/';'/';'/'};
    S4.HeatTransfer={num2str(Q_GycolCooling);num2str(Boiler);'/';'/';'/';'/';'/';'/'};
    S4.Temperature={'/';'/';'/';num2str(T_air);num2str(T_concentrated_milk);num2str(T_milk_powder);num2str(T_air_exit);num2str(T_water_exit)};
    S4.MassFlow={'/';'/';'/';num2str(m_air);num2str(m_concentrated_milk);num2str(m_milk_powder);num2str(m_air_sec);num2str(m_water)};
     S4.Cp={'/';'/';'/';num2str(cp_air_20);num2str(cp_concentrated_milk);num2str(cp_milk_powder);num2str(cp_air_20);num2str(cp_water_20)};
    
    DryerValues = struct2table(S4);

    writetable( DryerValues,'Spray_Dryer_Values_Germany.xls');
    Cost_Spray_dryer_Ge=double(Total);
end

end
disp 'Ending  <<Dryer section>> ';

%% Cleaning in place section

disp 'Starting  <<Cleaning section>> ';
for CleaningCountry=(1:2);
% 1) Mass&Energy balance :

T_water=25;

T_cip_e1=70;
m_cip_e1=4;
T_cip_e2=80;
T_cip_e3=75;
cp_water_70 = 4.191;
cp_water_25 = 4.180;
cp_water_75 = 4.194;
cp_water_80 = 4.198;

% If S-5 is isothermal, and cp==constant then :

T_recup_2=75;

syms m_water_var m_recup_2_var

[m_water m_recup_2]=solve(m_water_var*T_water + m_recup_2_var*T_recup_2 == m_cip_e1*T_cip_e1 , m_water_var + m_recup_2_var == m_cip_e1 , m_water_var, m_recup_2_var);

m_cip_e1=m_water+m_recup_2;
m_cip_e2=m_water+m_recup_2;
m_cip_e3=m_water+m_recup_2;
m_to_step_2=m_cip_e3-m_recup_2;

Q_STEAM=m_cip_e1*cp_water_75*(T_cip_e2-T_cip_e1);
Q_LAVAGE2=m_cip_e2*cp_water_75*(T_cip_e2-T_cip_e3);

% 2) Energy bill : 

ot=8000;
ep=[0.081 0.124];
gp=[0.038 0.05];

if CleaningCountry==1;

    Steam=(m_cip_e1*cp_water_75*(80-70))*ot*gp(1)/0.9;

    S5.Unit={'Mixer';'Steam';'LAVAGE2';'Water';'recup2';'cip-e1';'cip-e2';'cip-3';'to_step_2'};
    S5.Cost={'/';num2str(double(Steam));'/';'/';'/';'/';'/';'/';'/'};
    S5.HeatTransfer={'/';num2str(double(Q_STEAM));num2str(double(Q_LAVAGE2));'/';'/';'/';'/';'/';'/'};
    S5.Temperature={'/';'/';'/';num2str(T_water);num2str(T_recup_2);num2str(T_cip_e1);num2str(T_cip_e2);num2str(T_cip_e3);num2str(T_cip_e3)};
    S5.MassFlow={'/';'/';'/';num2str(double(m_water));num2str(double(m_recup_2));num2str(double(m_cip_e1));num2str(double(m_cip_e2));num2str(double(m_cip_e3));num2str(double(m_to_step_2))};
    
    CleaningValues = struct2table(S5);

    writetable( CleaningValues,'Cleaning_in_place_Values_France.xls');
    cost_cleaning_in_place_Fr=double(Steam);
end

if CleaningCountry==2;

    Steam=(m_cip_e1*cp_water_75*(80-70))*ot*gp(2)/0.9;

    S5.Unit={'Mixer';'Steam';'LAVAGE2';'Water';'recup2';'cip-e1';'cip-e2';'cip-3';'to_step_2'};
    S5.Cost={'/';num2str(double(Steam));'/';'/';'/';'/';'/';'/';'/'};
    S5.HeatTransfer={'/';num2str(double(Q_STEAM));num2str(double(Q_LAVAGE2));'/';'/';'/';'/';'/';'/'};
    S5.Temperature={'/';'/';'/';num2str(T_water);num2str(T_recup_2);num2str(T_cip_e1);num2str(T_cip_e2);num2str(T_cip_e3);num2str(T_cip_e3)};
    S5.MassFlow={'/';'/';'/';num2str(double(m_water));num2str(double(m_recup_2));num2str(double(m_cip_e1));num2str(double(m_cip_e2));num2str(double(m_cip_e3));num2str(double(m_to_step_2))};
    
    CleaningValues = struct2table(S5);

    writetable( CleaningValues,'Cleaning_in_place_Values_Germany.xls');
    cost_cleaning_in_place_Ge=double(Steam);

end

end
%% Cold&Heat Production 
disp 'Ending  <<Cleaning section>> ';
disp 'Starting  <<Storage section>> ';
for StorageCountry=(1:2);

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

Q_ColdStorage=500;
Q_E_31=m_hw_in*cp_hw_moy*(T_hw_out-T_hw_in);


COP=4;
ot=8000;
ep=[0.081 0.124];
gp=[0.038 0.05];

if StorageCountry==1;
    
   ColdStorage = 500/COP*ot*ep(1);

HotWater = m_hw_in*cp_hw_moy*(T_hw_out-T_hw_in)*ot*gp(1)/0.9;

Total = HotWater + ColdStorage; 
    
    S6.Unit={'ColdProduction';'HeatProduction';'Total'};
    S6.Cost={num2str(ColdStorage);num2str(HotWater);num2str(ColdStorage+HotWater)};
    S6.HeatTransfer={num2str(Q_ColdStorage);num2str(Q_E_31);'/'};
    
    StorageValues = struct2table(S6);

    writetable( StorageValues,'Heat&ColdProduction_France.xls');
    cost_cold_hot_Fr=ColdStorage+HotWater;
    
end


if StorageCountry==2;
    
   ColdStorage = 500/COP*ot*ep(2);

HotWater = m_hw_in*cp_hw_moy*(T_hw_out-T_hw_in)*ot*gp(2)/0.9;

Total = HotWater + ColdStorage; 
    
    S6.Unit={'ColdProduction';'HeatProduction';'Total'};
    S6.Cost={num2str(ColdStorage);num2str(HotWater);num2str(ColdStorage+HotWater)};
    S6.HeatTransfer={num2str(Q_ColdStorage);num2str(Q_E_31);'/'};
    
    StorageValues = struct2table(S6);

    writetable( StorageValues,'Heat&ColdProduction_Germany.xls');
cost_cold_hot_Ge=ColdStorage+HotWater;
end


end
disp 'Ending  <<Storage section>> ';

%% Evaporation section

%1 - Mass and energy balance

cp_water = 4.178; % trouvé dans la littérature
cp_milkPowd = 1.165;

T1 = 4;
T5 = 65;
T6 = 65;
T7 = 70;
T8 = 65;
T9 = 65;
T10 = 60;
T11 = 60;
T13 = 4;
T14 = 70;
T15 = 65;
T16 = 60;
T17 = 70;
T18 = 65;
T19 = 70;
T20 = 65;

DT_min = 5;
T12 = T1 + DT_min; %Dernier truc à vérifier
T3 = T16 - DT_min;
T4 = T15 - DT_min;
T_exit = T4 + DT_min;

m1 = 5;
m14 = 1.3;
m15 = 1.29;
m16 = 1.2;
m7 = m1-m14;
m9 = m7-m15;
m11 = m9-m16;

COP = 4;
efficiency = 0.9;
t = 8000; % Ecrit nul part

c_fuel = 0.038; % Pour la France
c_elec = 0.081; % Pour la France

cp_past = (0.121*m1/(m1))*cp_milkPowd + (1-(0.121*m1/(m1)))*cp_water;
cp_evap1 = (0.121*m1/(m1-m14))*cp_milkPowd + (1-(0.121*m1/(m1-m14)))*cp_water;
cp_evap2 = (0.121*m1/(m1-m14-m15))*cp_milkPowd + (1-(0.121*m1/(m1-m14-m15)))*cp_water;
cp_evap3 = (0.121*m1/(m1-m14-m15-m16))*cp_milkPowd + (1-(0.121*m1/(m1-m14-m15-m16)))*cp_water;

h_evap_200mbar = 2353;
h_evap_250mbar = 2348;
h_evap_310mbar = 2337;

m_evap_valve2 = (m7*cp_evap1*(T7-T8))/(h_evap_250mbar);
x_8 = m_evap_valve2/m7;
%fprintf('The vapor quality in 8 is equal to %g \n',x_8);
m_evap_valve3 = (m9*cp_evap2*(T9-T10))/(h_evap_200mbar);
x_10 = m_evap_valve3/m9;
%fprintf('The vapor quality in 10 is equal to %g \n \n',x_10);

Q_evap1 = m1*cp_past*(T14-T6) + m14*h_evap_310mbar;
Q_evap2 = m7*cp_evap1*(T15-T8) + (m15-m_evap_valve2)*h_evap_250mbar;
Q_evap3 = m9*cp_evap2*(T16-T10) + (m16-m_evap_valve3)*h_evap_200mbar;

m17 = Q_evap2/h_evap_310mbar;
m18 = Q_evap3/h_evap_250mbar;
m19 = m14 - m17;
m20 = m15 - m18;

% Pour les Q_HE, on a essayé de mettre le DT_min des deux côtés et on garde
% seulement la solution qui marche

Q_HE1 = m11*cp_evap3*(T11-T12);
T2 = T1 + Q_HE1/(m1*cp_past);

Q_HE2 = m1*cp_past*(T3-T2);

Q_HE3 = m1*cp_past*(T4-T3);

Q_HE4 = m19*h_evap_310mbar + m19*cp_water*(T19-T_exit);
T5 = T4 + Q_HE4/(m1*cp_past);

Q_boiler = Q_evap1/efficiency;
Q_ref = m11*cp_evap3*(T12-T13)/COP;

% 2 - Current energy bill

OC_boiler = Q_boiler*t*c_fuel;
OC_ref = Q_ref*t*c_elec;
OC_total = OC_boiler + OC_ref; %en euros !!!

TotalCost_Fr=Cost_Pasteurisation_Fr+Cost_Spray_dryer_Fr+cost_cleaning_in_place_Fr+cost_cold_hot_Fr
TotalCost_Ge=Cost_Pasteurisation_Ge+Cost_Spray_dryer_Ge+cost_cleaning_in_place_Ge+cost_cold_hot_Ge
TotalCost_Fr_Optimized=Cost_Pasteurisation_Fr_Optimized+Cost_Spray_dryer_Fr+cost_cleaning_in_place_Fr+cost_cold_hot_Fr
TotalCost_Ge_Optimized=Cost_Pasteurisation_Ge_Optimized+Cost_Spray_dryer_Ge+cost_cleaning_in_place_Ge+cost_cold_hot_Ge

disp 'Ended  <<Analysis>> ';