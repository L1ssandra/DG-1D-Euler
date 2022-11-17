function [Urhomod,Urhoumod,UEmod] = TVB_Limiter(Urho,Urhou,UE)

global Nx k bc gamma

% 边界条件
if bc == 1 % 周期
    UrhoR = [Urho(2:end,:);Urho(1,:)];
    UrhoL = [Urho(end,:);Urho(1:end - 1,:)];
    UrhouR = [Urhou(2:end,:);Urhou(1,:)];
    UrhouL = [Urhou(end,:);Urhou(1:end - 1,:)];
    UER = [UE(2:end,:);UE(1,:)];
    UEL = [UE(end,:);UE(1:end - 1,:)];
elseif bc == 2 % 进出口
    UrhoR = [Urho(2:end,:);Urho(end,:)];
    UrhoL = [Urho(1,:);Urho(1:end - 1,:)];
    UrhouR = [Urhou(2:end,:);Urhou(end,:)];
    UrhouL = [Urhou(1,:);Urhou(1:end - 1,:)];
    UER = [UE(2:end,:);UE(end,:)];
    UEL = [UE(1,:);UE(1:end - 1,:)];
end

Urhomod = zeros(Nx,k + 1);
Urhoumod = zeros(Nx,k + 1);
UEmod = zeros(Nx,k + 1);

for i = 1:Nx
    
    UR = [UrhoR(i,1) - Urho(i,1);UrhouR(i,1) - Urhou(i,1);UER(i,1) - UE(i,1)];
    UL = [Urho(i,1) - UrhoL(i,1);Urhou(i,1) - UrhouL(i,1);UE(i,1) - UEL(i,1)];
    Ux = [Urho(i,2);Urhou(i,2);UE(i,2)];
    
    u = Urhou(i,1)/Urho(i,1);
    p = (gamma - 1)*(UE(i,3) - 0.5*Urho(i,1)*u^2);
    H = (UE(i,3) + p)/Urho(i,1);
    c = sqrt(abs(gamma*p/Urho(i,1)));
    gamma1 = gamma - 1;
    
    R1 = [1;u - c;H - u*c];
    R2 = [1;u;0.5*u^2];
    R3 = [1;u + c;H + u*c];
     
    R = [R1,R2,R3] + 1e-10*eye(3);
    
    L = inv(R);
    
    Uxmod = R*minmod(L*Ux,L*UR,L*UL);
    
    Urhomod(i,1) = Urho(i,1);
    Urhomod(i,2) = Uxmod(1);
    
    Urhoumod(i,1) = Urhou(i,1);
    Urhoumod(i,2) = Uxmod(2);
    
    UEmod(i,1) = UE(i,1);
    UEmod(i,2) = Uxmod(3);
    
    if norm(Uxmod - Ux) == 0 && k > 1
        Urhomod(i,3:end) = Urho(i,3:end);
        Urhoumod(i,3:end) = Urhou(i,3:end);
        UEmod(i,3:end) = UE(i,3:end);
    end
    
end

end

function m = minmod(a,b,c)

global dim hx

M = 1; % TVB常数
m = a;

for i = 1:dim
    if a(i) > M*hx^2
        if sign(a(i)) == sign(b(i)) && sign(b(i)) == sign(c(i))
            m(i) = sign(a(i))*min([abs(a(i)),abs(b(i)),abs(c(i))]);
        else
            m(i) = 0;
        end
    end
end

end
    