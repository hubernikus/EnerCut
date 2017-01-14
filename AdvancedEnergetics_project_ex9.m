% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
%
% Projet - Advanced ENergetics
%
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
clc; close all; clear variables;
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
%% Task 1
% Constants
T_amb = 20; %[C] Ambient temperature
COP = 4;

% Define Dataset
N_flows = 8;

Tin = [75 65 10 25 20 0 10 90];
Tin = Tin(1:N_flows);   

Tout = [25 25 3 80 70 0.1 15 60];
Tout = Tout(1:N_flows);

% Enthalpy
Q = [665 588 300 836 924 0 0 0];
Q = Q(1:N_flows);

% Enthalpy hot utitlity
deltaT = [5 5 2 5 5 1 2 5];
deltaT = deltaT(1:N_flows);

alpha = [0.5 0.5 1.5 0.5 0.5 2 2 1.5];
alpha = alpha(1:N_flows);

N = length(Tin);

hotFlow = (Tin-Tout) > 0;

% T_star_hot = zeros(N,2);
Tin_dt = zeros(1,N);
Tout_dt = zeros(1,N);

for i = 1:length(Tin)
    if hotFlow(i)
        Tin_dt(i) = Tin(i) - deltaT(i);
        Tout_dt(i) = Tout(i) - deltaT(i);
    else    % cold flow
        Tin_dt(i) = Tin(i) + deltaT(i);
        Tout_dt(i) = Tout(i) + deltaT(i);
    end
end

% Sort the temperature 
[T_sort,~,ind_T] = unique([Tin_dt,Tout_dt]);
T_sort = flip(T_sort); ind_T = length(T_sort)+1- ind_T'; % Reverse order (descending)
N_sortT= length(T_sort);

% Enthalpy - Temperature Diagramm
hotCurve = zeros(1, N_sortT);
coldCurve = zeros(1, N_sortT);

for i = 1:5 % The number of streams without utilities & refrigiration cycle
    if hotFlow(i)
        ind_top = ind_T(i+N);
        T_top = T_sort(ind_top);
        T_bot = T_sort(ind_T(i));
        for jj = 1:ind_top
            hotCurve(jj) = hotCurve(jj) + Q(i)*min(1,(T_sort(jj)-T_top)/(T_bot-T_top));
        end
    else
        ind_top = ind_T(i);
        T_top = T_sort(ind_top);
        T_bot = T_sort(ind_T(i+N));
        for jj = 1:ind_top
            coldCurve(jj) = coldCurve(jj) + Q(i)*min(1,(T_sort(jj)-T_top)/(T_bot-T_top));
        end
    end
end

% Allign hot and cold stream so that they touch at 10 deg.
ind10 = find(T_sort == 12);
coldCurve = coldCurve + hotCurve(10)-coldCurve(10);

% Define missing Enthalpies
ind20 = find(T_sort == 20); % Find 20 deg to exchange heat.

ind = 6; % Vecotr index of referigiation cycle
Q(ind) = coldCurve(end)-hotCurve(end);
ind_top = ind_T(ind); T_top = T_sort(ind_top); T_bot = T_sort(ind_T(ind+N));
for jj = 1:ind_top
    coldCurve(jj) = coldCurve(jj) + Q(ind)*min(1,(T_sort(jj)-T_top)/(T_bot-T_top));
end
coldCurve = coldCurve + hotCurve(10)-coldCurve(10);

ind = 7; % Vector index of hot utility
Q(ind) = max(hotCurve(1:ind10)-coldCurve(1:ind10)); % Hot Utility
ind_top = ind_T(ind); T_top = T_sort(ind_top); T_bot = T_sort(ind_T(ind+N));
for jj = 1:ind_top
    coldCurve(jj) = coldCurve(jj) + Q(ind)*min(1,(T_sort(jj)-T_top)/(T_bot-T_top));    
end

% Cold utility
ind = 8;
Q(ind) = coldCurve(1)-hotCurve(1);
ind_top = ind_T(ind+N); T_top = T_sort(ind_top); T_bot = T_sort(ind_T(ind));
for jj = 1:ind_top
    hotCurve(jj) = hotCurve(jj) + Q(ind)*min(1,(T_sort(jj)-T_top)/(T_bot-T_top));
end

dT_sort = T_sort(1:end-1) - T_sort(2:end);

TempTable = table(coldCurve', hotCurve', 'VariableNames', ... 
            {'coldCurve','hotCurve'});


% % Plot Data 
figure;
plot(hotCurve, T_sort, 'r'); hold on;
plot(coldCurve, T_sort, 'b'); hold on;
legend('Hot streams', 'Cold streams', 'Location', 'east');
xlabel('Heat load [kW]'); ylabel('Temeperature [C]');
xlim([0,hotCurve(1)]); ylim([0,1000]);
title('Composite Curve');
grid on; ylim([0,100])
print('fig_copositeCurve','-dpng')

%% Task 2
% Calulate 
E_T = abs(Q./(Tin-Tout));

M_cp = zeros(1,N_sortT-1);

for i = 1:N
    Mi = zeros(1,N_sortT-1);
    if hotFlow(i)
        Mi(ind_T(i):ind_T(i+N)-1) = E_T(i);
    else
        Mi(ind_T(i+N):ind_T(i)-1) = -1*E_T(i);
    end
    M_cp = M_cp + Mi;
end

Q_r = dT_sort.*M_cp;

R_r0 = zeros(1,N_sortT);    

for i = 2:N_sortT
   R_r0(i) = R_r0(i-1)+Q_r(i-1);
end

R_r = R_r0 - min(R_r0);

Mcp = Q./abs(Tin-Tout);

Q_above = zeros(1,N);
Q_below = zeros(1,N);

T_pinch = 30;
for i = 1:N
    if hotFlow(i)
        Q_above(i) = min(max((Tin_dt(i)-T_pinch)*Mcp(i),0),Q(i));
        Q_below(i) = Q(i)-Q_above(i);
    else
        Q_above(i) = min(max((Tout_dt(i)-T_pinch)*Mcp(i), 0),Q(i));
        Q_below(i) = Q(i)-Q_above(i);
    end
end

Hu2 = sum(Q_above(4:5))-sum(Q_above(1:3));
Cu2 = sum(Q_below(1:3))- sum(Q_below(4:5));
% Generate Table
TemperatureTable = table(Tin_dt', Tout_dt', Q',Mcp', Q_above', Q_below', ...
                    'VariableNames', {'T_in_dt','T_out_dt', 'Q','M_cp','Q_above','Q_below'}, ...
                    'RowNames', {'H1','H2','H3','C1','C2','Re','Cw','Hu'})

NetWork_Table = table(T_sort', [0 dT_sort]', [0 M_cp]', [0 Q_r]', R_r0', R_r', ...
                'VariableNames', {'TC','deltaT','M_cp','Q_r','R_r0','R_r'})

figure;
plot(R_r, T_sort,'r'); hold on;
plot([R_r(1) R_r(1)],[0,T_sort(1)],'--k'); hold on;
plot([R_r(end) R_r(end)],[0 T_sort(end)],'--k'); 
xlabel('Heat Load [kW]'); ylabel('Temperature [C]');
legend('All streams'); grid on;
title('Grand Composite Curve');
print('fig_grandCompositeCurve','-dpng')
            

%% Calculation Cost
price_fuel = [0.050, 0.038];      % [EUR / kWh]   -> price Germany: 0.50  // price France: 0.038
price_elec = [0.124, 0.081];      % [EUR / kWh]   -> price Germany: 0.124  // price France: 0.081

i = 2; % 1-Germany / 2-France
price = [price_fuel(i), price_elec(i)];


price_wat = 0.01; % [0.01 CHF]
% Q_R is rectified to see in the end what energy amount has to added using
% CU & HU

% Define System1
hot = [1,2,2,1];
cold = [1,2,2,2];
abovePinch = [1,1,0,0];

Q_R1 = [Q_above; Q_below]; 
Mcp_R1 = Mcp;

% Run calcluations
for i = 1:length(cold)
    [Q_tran1(i), A1(i), CostUnit1(i),  investmentCost1(i), Q_R1, Mcp_R1, ~] = ...
                        costCalcul(hot(i), cold(i), abovePinch(i), Q_R1, ...
                        T_pinch, Tin_dt, deltaT, alpha, Mcp_R1);
end

Q_HU1 = sum(Q_R1(1,4:5));
Q_CU1 = sum(Q_R1(2,1:3));
if (round(Q_HU1 -Q(8)))
    warning('Hot utility: Q_8=%d, Q_H=%d',Q(8),Q_HU1)
elseif(round(Q_CU1 - sum(Q(6:7))))
    warning('Cold utility: Q_67=%d, Q_CU=%d',sum(Q(6:7)),Q_CU1)
else
    fprintf('Cold utiility: CU & HU all good \n')
end

%NetWork_Table 
NetWork_Table = table(Q_tran1', A1', CostUnit1', investmentCost1', ...
                'VariableNames', {'Q_i','A','CostUnit', 'InvestmentCost'})

% Define System2
hot2 = [2,1,1,2];
cold2 = [1,2,2,2];
abovePinch2 = [1,1,0,0];

Q_R2 = [Q_above; Q_below]; 
Mcp_R2 = Mcp;

% Run calcluations
for i = 1:length(cold)
    [Q_tran2(i),A2(i), CostUnit2(i),  investmentCost2(i), Q_R2, Mcp_R2, ~] = ...
                        costCalcul(hot2(i), cold2(i), abovePinch2(i), Q_R2, ...
                        T_pinch, Tin_dt, deltaT, alpha, Mcp_R2);
end

Q_HU2 = sum(Q_R2(1,4:5));
Q_CU2 = sum(Q_R2(2,1:3));
if (round(Q_HU2 -Q(8)))
    warning('Hot utility: Q_8=%d, Q_H=%d',Q(8),Q_HU2)
elseif(round(Q_CU2 - sum(Q(6:7))))
    warning('Cold utility: Q_67=%d, Q_CU=%d',sum(Q(6:7)),Q_CU1)
else
    fprintf('Cold utiility: CU & HU all good \n')
end

%NetWork_Table 
NetWork_Table = table(A2', CostUnit2', investmentCost2', ...
                'VariableNames', {'A','CostUnit','InvestmentCost'});
            

% Print to Latex
fileID1 = fopen('HEN_opt1.tex','w');
formatSpec = 'E%d & %3.0f & %3.0f & %6.0f & %6.0f \\\\ \n';
fprintf(fileID1, formatSpec, [(1:length(Q_tran1));round(Q_tran1); round(A1);...
                round(CostUnit1); round(investmentCost1)] );
fclose(fileID1);            

% Print to Latex
fileID2 = fopen('HEN_opt2.tex','w');
formatSpec = 'E%d & %3.0f & %3.0f & %6.0f & %6.0f \\\\ \n';
fprintf(fileID2, formatSpec, [(1:length(Q_tran2));round(Q_tran2); round(A2);...
                round(CostUnit2); round(investmentCost2)] );
fclose(fileID2);    


% Final Option
% Define System1
hot = [1,2,1,2];
cold = [1,2,2,2];
abovePinch = [1,1,0,0];

Q_R = [Q_above; Q_below]; 
Mcp_R = Mcp;

% Run calcluations
for i = 1:length(cold)
    [Q_tran(i), A(i), CostUnit(i),  investmentCost(i), Q_R, Mcp_R, ~] = ...
                        costCalcul(hot(i), cold(i), abovePinch(i), Q_R, ...
                        T_pinch, Tin_dt, deltaT, alpha, Mcp_R);
end

Q_HU = sum(Q_R(1,4:5));
Q_CU = sum(Q_R(2,1:3));
if (round(Q_HU -Q(8)))
    warning('Hot utility: Q_8=%d, Q_H=%d',Q(8),Q_HU)
elseif(round(Q_CU - sum(Q(6:7))))
    warning('Cold utility: Q_67=%d, Q_CU=%d',sum(Q(6:7)),Q_CU)
else
    fprintf('Cold utiility: CU & HU all good \n')
end

fileID = fopen(['costAnnual.tex'],'w');
fprintf(fileID,'\\begin{tabular}{l r r r r} \n \\hline \n');
fprintf(fileID, ' & Investment & CAPEX & OPEX & Total \\\\ \n');
fprintf(fileID, ' & [CHF] & [CHF/year] & [CHF/year] & [CHF/year] \\\\  \\hline \n'); 

Q_ref = Q(3); % Heat refigiration cycle
annualCost(0, sum(Q(4:5)),sum(Q(1:2)),Q_ref, fileID,'No Exchanger', price);
annualCost(investmentCost, sum(Q_HU),sum(Q_CU),Q_ref, fileID,'Option A-B', price);

hot = 1; cold = 2; abovePinch = true; Q_max = Q_tran(3);
[Q_transf1,A_2,~,invCost1,~,~,T_out2] = costCalcul(hot, cold, abovePinch, Q, T_pinch, ...
                                Tin_dt, deltaT, alpha, Mcp, Q_max)
                            
hot = 2; cold = 2; abovePinch = true; Q_max = Q_tran(2);
[Q_transf2,~,~,invCost2,~,~,~] = costCalcul(hot, cold, abovePinch, Q, T_out2(2), ...
                                Tin_dt, deltaT, alpha, Mcp, Q_max)

investmentCost_new = investmentCost(1) +invCost1 +invCost2;                            
Q_HU_new = sum(Q(4:5))- Q_tran(1)-Q_transf1-Q_transf2;
Q_CU_new = sum(Q(1:2))- Q_tran(1)-Q_transf1-Q_transf2;
annualCost(investmentCost_new, Q_HU_new, Q_CU_new, Q_ref, fileID,'3 Exchanger', price);

hot = 2; cold = 2; abovePinch = true; Q_max = Q(2)*(Tin_dt(2)-Tin_dt(2+3))/(Tin_dt(2)-Tout_dt(2))
[Q_transf1,~,~,invCost1,~,~,T_out2] = costCalcul(hot, cold, abovePinch, Q, Tin(2+3), ...
                                Tin_dt, deltaT, alpha, Mcp, Q_max)
                            
investmentCost_new = investmentCost(1) +invCost1;
Q_HU_new = sum(Q(4:5))- Q_tran(1)-Q_transf1;
Q_CU_new = sum(Q(1:2))- Q_tran(1)-Q_transf1;
annualCost(investmentCost_new, Q_HU_new, Q_CU_new, Q_ref, fileID,'2 Exchanger', price);

fprintf(fileID, '\\hline \n \\end{tabular} \n');
fclose(fileID);

% opCost_ref = W_ref*price_elec*opTime

%% Task 3
% Pressure Ratio Evaporator
r = 4;                          % p_cond/p_eva
eff_compr = 0.82;                % h2 = (h2s-1)*eff +h1
T_evap = Tin(5)
p_evap0_R134A = 298200;          % [Pa] 
p_cond_R134A = r*p_evap0_R134A % [Pa]

T_cond = 46.04 % C -> calculated with pressure

h_evap2 = 250.46; %[kJ/kg]
s_evap2 = 0.9315; %[kJ/kg K]
    
h_cond1 = 307.58; %[kJ/kg]
ds_cond1 = 1.013; %[ kJ/kg K]

h_cond2 = 112.77; %[kJ/kg]
s_cond2 = 0.409; %[kJ/kg K]

h_evap1 = h_cond2; % [kJ/kg]
s_evap1 = 0.4269; % [kJ/kg]

% Exergy - Balance
Q_evap = Q(5)

% Crossflow
dt_Out = Tin(3)-Tout(6); % cold
dt_In = Tout(3) - Tin(6); % hot
lmdt_evap = (dt_Out-dt_In)/log(dt_Out/dt_In);
m_ref = Q_evap/(h_evap2-h_evap1);
A_evap = Q_evap/(alpha(5)*lmdt_evap);

% cross flow
dt_Out = T_cond - T_amb; %
dt_In = T_cond - (T_amb + 5); % Air is heated by 5 degrees
lmdt_cond = (dt_Out-dt_In)/log(dt_Out/dt_In);

Q_cond = m_ref*(h_cond1-h_cond2)
A_cond = Q_cond/(alpha(5)*lmdt_cond);

W_ref = Q_cond - Q_evap

% Exergy analysis
X_ref = Q_evap*T_amb/T_evap-Q_cond*T_amb/T_cond + W_ref;

