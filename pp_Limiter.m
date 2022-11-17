function [Urho,Urhou,UE] = pp_Limiter(Urho,Urhou,UE)

epsilon = 1e-13;

global Nx phiValueR phiValueL phiValue0

% 限制 rho
for i = 1:Nx
    % 计算三个积分点上的值
    rhoR = Urho(i,1)*phiValueR(1) + Urho(i,2)*phiValueR(2) + Urho(i,3)*phiValueR(3);
    rho0 = Urho(i,1)*phiValue0(1) + Urho(i,2)*phiValue0(2) + Urho(i,3)*phiValue0(3);
    rhoL = Urho(i,1)*phiValueL(1) + Urho(i,2)*phiValueL(2) + Urho(i,3)*phiValueL(3);
    
    rhomin = min([rhoR,rho0,rhoL]);
    
    theta = min( (Urho(i,1) - epsilon)/(Urho(i,1) - rhomin),1 );
    
    Urho(i,2:end) = theta*Urho(i,2:end);
end

% 限制p
for i = 1:Nx
    
    Ubar = [Urho(i,1);Urhou(i,1);UE(i,1)];
    
    rhoL = Urho(i,1)*phiValueL(1) + Urho(i,2)*phiValueL(2) + Urho(i,3)*phiValueL(3);
    rhouL = Urhou(i,1)*phiValueL(1) + Urhou(i,2)*phiValueL(2) + Urhou(i,3)*phiValueL(3);
    EL = UE(i,1)*phiValueL(1) + UE(i,2)*phiValueL(2) + UE(i,3)*phiValueL(3);
    
    U1 = [rhoL;rhouL;EL];
    t1 = bisectP(U1,Ubar);
    s1 = (1 - t1)*Ubar + t1*U1;
    S1 = norm(s1 - Ubar)/(norm(U1 - Ubar) + epsilon);
    
    rho0 = Urho(i,1)*phiValue0(1) + Urho(i,2)*phiValue0(2) + Urho(i,3)*phiValue0(3);
    rhou0 = Urhou(i,1)*phiValue0(1) + Urhou(i,2)*phiValue0(2) + Urhou(i,3)*phiValue0(3);
    E0 = UE(i,1)*phiValue0(1) + UE(i,2)*phiValue0(2) + UE(i,3)*phiValue0(3);
    
    U2 = [rho0;rhou0;E0];
    t2 = bisectP(U2,Ubar);
    s2 = (1 - t2)*Ubar + t2*U2;
    S2 = norm(s2 - Ubar)/(norm(U2 - Ubar) + epsilon);
    
    rhoR = Urho(i,1)*phiValueR(1) + Urho(i,2)*phiValueR(2) + Urho(i,3)*phiValueR(3);
    rhouR = Urhou(i,1)*phiValueR(1) + Urhou(i,2)*phiValueR(2) + Urhou(i,3)*phiValueR(3);
    ER = UE(i,1)*phiValueR(1) + UE(i,2)*phiValueR(2) + UE(i,3)*phiValueR(3);
    
    U3 = [rhoR;rhouR;ER];
    t3 = bisectP(U3,Ubar);
    s3 = (1 - t3)*Ubar + t3*U3;
    S3 = norm(s3 - Ubar)/(norm(U3 - Ubar) + epsilon);
    
    theta = min([S1,S2,S3]);
    
    Urho(i,2:end) = theta*Urho(i,2:end);
    Urhou(i,2:end) = theta*Urhou(i,2:end);
    UE(i,2:end) = theta*UE(i,2:end);
    
end

end

function t = bisectP(U,Ubar)

epsilon = 1e-13;
t0 = 0; t1 = 1;

if pressure(U) > epsilon
    t = 1;
else
    error = abs(pressure(Ubar) - epsilon);
    N = 0;
    while error > 1e-14 && N < 50
        t = 0.5*(t0 + t1);
        Ut = (1 - t)*Ubar + t*U;
        if pressure(Ut) > epsilon
            t0 = t;
        else
            t1 = t;
        end
        N = N + 1;
    end
    t = t0;
end

end