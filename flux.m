function Fhat = flux(FL,FR,rhoL,rhouL,EL,rhoR,rhouR,ER,SR,SL)

global gamma type

UR = [rhoR;rhouR;ER]; UL = [rhoL;rhouL;EL];

switch type
    
    case 1
       %% L-F Flux
        Fhat = 0.5*(FR + FL - max(abs(SR),abs(SL))*(UR - UL));
        
    case 2
       %% HLL Flux
        if SL >= 0
            Fhat = FL;
        elseif SR <= 0
            Fhat = FR;
        else
            Fhat = (SR*FL - SL*FR + SL*SR*(UR - UL))/(SR - SL);
        end
        
    case 3
       %% HLLC Flux
        if SL >= 0
            Fhat = FL;
        elseif SR <= 0
            Fhat = FR;
        else
            
            uL = rhouL/rhoL;
            uR = rhouR/rhoR;
            
            pL = (gamma - 1)*(EL - 0.5*rhoL*uL^2);
            pR = (gamma - 1)*(ER - 0.5*rhoR*uR^2);
            
            Sstar = (pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR))/(rhoL*(SL - uL) - rhoR*(SR - uR));
            
            Dstar = [0;1;Sstar];
            
            if SL < 0 && Sstar >= 0
                Fhat = ( Sstar*(SL*UL - FL) + SL*(pL + rhoL*(SL - uL)*(Sstar - uL))*Dstar )/(SL - Sstar);
            else
                Fhat = ( Sstar*(SR*UR - FR) + SR*(pR + rhoR*(SR - uR)*(Sstar - uR))*Dstar )/(SR - Sstar);
            end
        end
end

end