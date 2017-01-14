function [] =  annualCost(investmentCost, Q_HU, Q_CU, Q_ref, fileID, modification, price)


price_fuel = price(1);      % [EUR / kWh] -> price Germany: 0.050  // price France: 0.038
price_elec = price(2);      % [EUR / kWh] -> price Germany: 0.124  // price France: 0.081

price_wat = 0.01; % [CHF/l]    

%Operating time per year
opTime = 8000;      % [h/year]

% Efficiency hot heating
effHeat = 0.85; 

% Hot Utility
OPEX_HU = Q_HU*opTime*price_fuel/effHeat;

% Cold Utility
cp_wat = 4.178;                                  % [kg/kg K]    
rho_wat = 996;                                  % [kg/m^3]

Tin_wat = 10;   % [C]
Tout_wat = 15;  % [C] 

M_CU = Q_CU/(cp_wat*(Tout_wat-Tin_wat));
OPEX_CU = M_CU*3600/rho_wat*price_wat*opTime; % 3600 -> seconds per hour]

% Efficiency Heatpump
COP = 4;

% Refigiration Cycle
OPEX_RE = Q_ref*opTime*price_elec/COP;

OPEX = OPEX_HU + OPEX_CU + OPEX_RE;

% Amortization Time
tAmort = 8.55; % [years]

CAPEX = sum(investmentCost/tAmort);

fprintf(fileID, '%s & %d & %d & %d & %d \\\\ \n', ...
            modification, round(sum(investmentCost)), round(CAPEX), ...
            round(OPEX), round(CAPEX+OPEX)); 
end