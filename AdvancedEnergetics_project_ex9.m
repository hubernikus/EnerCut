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
            {'coldCurve','hotCurve'})

% Plot Data 
% figure;
% plot(hotCurve, T_sort, 'r'); hold on;
% plot(coldCurve, T_sort, 'b'); hold on;
% legend('Hot streams', 'Cold streams', 'Location', 'east');
% xlabel('Heat load [kW]'); ylabel('Temeperature [C]');
% xlim([0,hotCurve(1)]); ylim([0,1000]);
% title('Composite Curve');
% grid on; ylim([0,100])
% print('fig_copositeCurve','-dpng')

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

Hu2 = sum(Q_above(4:5))-sum(Q_above(1:3))
Cu2 = sum(Q_below(1:3))- sum(Q_below(4:5))
% Generate Table
TemperatureTable = table(Tin_dt', Tout_dt', Q',Mcp', Q_above', Q_below', ...
                    'VariableNames', {'T_in_dt','T_out_dt', 'Q','M_cp','Q_above','Q_below'}, ...
                    'RowNames', {'H1','H2','H3','C1','C2','Re','Cw','Hu'})

NetWork_Table = table(T_sort', [0 dT_sort]', [0 M_cp]', [0 Q_r]', R_r0', R_r', ...
                'VariableNames', {'TC','deltaT','M_cp','Q_r','R_r0','R_r'})

% figure;
% plot(R_r, T_sort,'r'); hold on;
% plot([R_r(1) R_r(1)],[0,T_sort(1)],'--k'); hold on;
% plot([R_r(end) R_r(end)],[0 T_sort(end)],'--k'); 
% xlabel('Heat Load [kW]'); ylabel('Temperature [C]');
% legend('All streams'); grid on;
% title('Grand Composite Curve');
% print('fig_grandCompositeCurve','-dpng')
            

%% Calculation Cost
price_fuel = [0.50, 0.038];      % [EUR / kWh]   -> price Germany: 0.50  // price France: 0.038
price_elec = [0.124, 0.081];      % [EUR / kWh]   -> price Germany: 0.124  // price France: 0.081

price_fuel = price_fuel(1);
price_elec = price_elec(1);


price_wat = 0.01; % [0.01 CHF]
% Q_R is rectified to see in the end what energy amount has to added using
% CU & HU

% Define System1
hot = [1,2,2,1];
cold = [4,5,5,5];
abovePinch = [1,1,0,0];

Q_R = [Q_above; Q_below]; 
Mcp_R = Mcp;

% Run calcluations
for i = 1:length(cold)
    disp(i)
    [Q_tran(i), A(i), CostUnit(i),  investmentCost(i), Q_R, Mcp_R] = ...
                        costCalcul(hot(i), cold(i), abovePinch(i), Q_R, ...
                        T_pinch, Tin, Tout, deltaT, alpha, Mcp_R);
end

Q_HU = sum(Q_R(1,4:5));
Q_CU = sum(Q_R(2,1:3));
if (round(Q_HU -Q(8)))
    warning('Hot utility: Q_8=%d, Q_H=%d',Q(8),Q_HU)
elseif(round(Q_CU - sum(Q(6:7))))
    warning('Cold utility: Q_67=%d, Q_CU=%d',sum(Q(6:7)),Q_CU)
else
    fprintf('Cold utiility: CU & HU all good')
end

%NetWork_Table 
NetWork_Table = table(Q_tran', A', CostUnit', investmentCost', ...
                'VariableNames', {'Q_i','A','CostUnit','InvestmentCost'})

% Print to Latex
fileID1 = fopen('HEN_opt1.tex','w');
formatSpec = 'E%d & %3.0f & %3.0f & %6.0f & %6.0f \\\\ \n';
fprintf(fileID1, formatSpec, [(1:length(Q_tran));round(Q_tran); round(A,0);...
                round(CostUnit,-3); round(investmentCost, -3)] );
fclose(fileID1);            
            
% Define System2
hot = [2,1,1,2];
cold = [4,5,5,5];
abovePinch = [1,1,0,0];

Q_R = [Q_above; Q_below]; 
Mcp_R = Mcp;
% Run calcluations
for i = 1:length(cold)
    [Q_tran(i),A(i), CostUnit(i),  investmentCost(i), Q_R, Mcp_R] = ...
                        costCalcul(hot(i), cold(i), abovePinch(i), Q_R, ...
                        T_pinch, Tin, Tout, deltaT, alpha, Mcp_R);
end

Q_HU = sum(Q_R(1,4:5));
Q_CU = sum(Q_R(2,1:3));
if (round(Q_HU -Q(8)))
    warning('Hot utility: Q_8=%d, Q_H=%d',Q(8),Q_HU)
elseif(round(Q_CU - sum(Q(6:7))))
    warning('Cold utility: Q_67=%d, Q_CU=%d',sum(Q(6:7)),Q_CU)
else
    fprintf('Cold utiility: CU & HU all good')
end

% Print to Latex
fileID2 = fopen('HEN_opt2.tex','w');
formatSpec = 'E%d & %3.0f & %3.0f & %6.0f & %6.0f \\\\ \n';
fprintf(fileID2, formatSpec, [(1:length(Q_tran));round(Q_tran); round(A,0);...
                round(CostUnit,-3); round(investmentCost, -3)] );
fclose(fileID2);            

%NetWork_Table 
NetWork_Table = table(A', CostUnit', investmentCost', ...
                'VariableNames', {'A','CostUnit','InvestmentCost'})
            

%Operating time per year
opTime = 8000; % h/year 

%Hot Utility
eff_heat = 0.85;

Q_HU = Q(8); %245+145;
Q_gas = Q_HU*eff_heat./price_fuel;

opCost_HU = Q_gas*price_fuel*opTime;

% Cold Utility
Q_CU = Q(7);
cp_wat =4.178;  
rho_wat = 996;              %[kg/m^3]
M_CU = Q_CU/(cp_wat*(Tout(7)-Tin(7)));
opCost_CU = M_CU*3600/rho_wat*price_wat*opTime; % 3600s/h

Q_ref = Q(6);

W_ref = Q_ref/4;

% opCost_ref = W_ref*price_elec*opTime



%% Task 3
% Pressure Ratio Evaporator
r = 4;                          % p_cond/p_eva
eff_compr = 0.82;                % h2 = (h2s-1)*eff +h1
T_evap = Tin(5);
p_evap0_R134A = 298200;          % [Pa] 
p_cond_R134A = r*p_evap0_R134A; % [Pa]

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
lmdt_evap = (dt_Out-dt_In)/log(dt_Out/dt_In)
m_ref = Q_evap/(h_evap2-h_evap1);
A_evap = Q_evap/(alpha(5)*lmdt_evap)

% cross flow
dt_Out = T_cond - T_amb; %
dt_In = T_cond - (T_amb + 5); % Air is heated by 5 degrees
lmdt_cond = (dt_Out-dt_In)/log(dt_Out/dt_In)

Q_cond = m_ref*(h_cond1-h_cond2)
A_cond = Q_cond/(alpha(5)*lmdt_cond)

W_ref = Q_cond - Q_evap;

% Exergy analysis
X_ref = Q_evap*T_amb/T_evap-Q_cond*T_amb/T_cond + W_ref


%% Functions
function [Q_ex, A, CostUnit,  investmentCost, Q_R, Mcp] ... 
    = costCalcul(hot, cold, abovePinch, Q_R, T_pinch, Tin, Tout, deltaT, alpha, Mcp)
Q_aboveR = Q_R(1,:);
    Q_belowR = Q_R(2,:);

    if(abovePinch) 
        Q_ex = min(Q_aboveR(hot), Q_aboveR(cold));

        Q_aboveR(hot) = Q_aboveR(hot)-Q_ex;
        Q_aboveR(cold) = Q_aboveR(cold)-Q_ex;

        dt1 = Tin(hot)-(Q_ex/Mcp(cold)+T_pinch); % For this example we start at pinch
        dt2 = max(Tin(hot)- Q_ex/Mcp(hot), T_pinch) - (T_pinch);
    else
        Mcp_temp = 0;
        Q_ex = min(Q_belowR(hot), Q_belowR(cold));
        if(Mcp(hot) < Mcp(cold)) % Split stream
            Q_ex = Q_ex*Mcp(hot)/Mcp(cold);
            Mcp_temp = Mcp(cold);
            Mcp(cold) = Mcp(hot);
            
        end
        Q_belowR(hot) = Q_belowR(hot)-Q_ex; % For this example we start at pinch
        Q_belowR(cold) = Q_belowR(cold)-Q_ex;

        dt1 = (T_pinch-(Q_ex/Mcp(cold)+Tin(cold)));
        dt2 = (T_pinch -Q_ex/Mcp(hot)) - Tin(cold);
        if(Mcp_temp) % Assign resting Mcp
            Mcp(cold) = Mcp_temp -Mcp(hot);
        end
    end
    Q_R = [Q_aboveR; Q_belowR]; 
    
    dtMin = deltaT(hot)+deltaT(cold);
    if(or(dt2< dtMin, dt1<dtMin))
       warning('Exchange hardly possible. dt2 too small.')
       fprintf('dt1=%d, dt2=%d ;; dtMin=%d\n',dt1,dt2,dtMin)
    end;
    
    % Logarithmic mean temperature difference (LMTD) in correcteddomain
    if dt1 == dt2
        lmtd = dt1;
    else
        lmtd = (dt1-dt2)/log(dt1/dt2);
    end
    
    U = 1/(1/alpha(hot)+1/alpha(cold)); % Global heat transfer coefficient
    A = Q_ex/(lmtd*U); % Heat exchanger area

    k1 = 3.8528; k2 = 0.4242; % Reference year 1998
    I_t = 585.7; % most recent available 2011
    I_ref = 389.5; % 1998
    CostUnit = I_t/I_ref*10^(k1+k2*log(A));% Investment cost
    
    nYears = 20; % Lifetime of the installation
    intRate = 0.08; % Interest                                                                                                                                                                                                                                                                                                                      
    exch_chfUSD = 0.99; % @ 13th Dec 2016
    F_BM = 4.74;% Bare module factor

    purchaseCost = F_BM*sum(CostUnit)*exch_chfUSD;

    investmentCost = purchaseCost * intRate*(1+intRate)^nYears/((1+intRate)^nYears -1);
end