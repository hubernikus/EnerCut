function [Q_ex, A, CostUnit,  investmentCost, Q_R, Mcp, Tout] ... 
    = costCalcul(hot, cold, abovePinch, Q_R, T_limit, Tin_dt, deltaT, alpha, Mcp, Q_max)
% COSTCALCUL Caluclation of cost for the heat exchanger network
cold = cold +3;
    Q_aboveR = Q_R(1,:);
    if size(Q_R,1) == 1 % no Q_below assigned -> use twice the same
        Q_belowR = Q_aboveR;
    else
        Q_belowR = Q_R(2,:);
    end
    
    if(abovePinch) 
        if nargin == 10 % Q_max is for HEN optimisation and limits Q 
            Q_ex = Q_max;
        else
            Q_ex = min(Q_aboveR(hot), Q_aboveR(cold));
        end
        
        Q_aboveR(hot) = Q_aboveR(hot)-Q_ex;
        Q_aboveR(cold) = Q_aboveR(cold)-Q_ex;

        Tout(1) = Tin_dt(hot)- Q_ex/Mcp(hot);
        Tout(2) = Tin_dt(cold) + Q_ex/Mcp(cold);
        
        dt1 = Tin_dt(hot)- Tout(2); % For this example we start at pinch
        dt2 = Tout(1) - T_limit;
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

        Tout(1) = Tin_dt(hot)- Q_ex/Mcp(hot);
        Tout(2) = Tin_dt(cold) + Q_ex/Mcp(cold);
        
        dt1 = Tin_dt(hot)- Tout(2); % For this example we start at pinch
        dt2 = Tout(1) - T_limit;
        
        %dt1 = (T_limit-(Q_ex/Mcp(cold)+Tin(cold)));
        %dt2 = (T_limit -Q_ex/Mcp(hot)) - Tin(cold);
        if(Mcp_temp) % Assign resting Mcp
            Mcp(cold) = Mcp_temp -Mcp(hot);
        end
    end
    Q_R = [Q_aboveR; Q_belowR]; 
    
    dtMin = deltaT(hot)+deltaT(cold);
    if(or(dt2< 0, dt1<0))
       warning('Exchange hardly possible. dt too small.')
       fprintf('dt1=%d, dt2=%d, dtMin=%d\n',dt1,dt2,dtMin)
    end;
    
    % Logarithmic mean temperature difference (LMTD) in corrected domain
    if dt1 == dt2
        lmtd = dt1;
    else
        lmtd = (dt1-dt2)/log((dt1+dtMin)/(dt2+dtMin));
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