function [dUrho,dUrhou,dUE] = Lh(Urho,Urhou,UE)

global hx1 Nx k quad phiValueInt phiValueIntx phiValueL phiValueR weight mm bc

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

dUrho = zeros(Nx,k + 1);
dUrhou = zeros(Nx,k + 1);
dUE = zeros(Nx,k + 1);


for i = 1:Nx
    % 体积分 \int_Ii f(u)*v_x dx 
    Qrho = Urho(i,1)*phiValueInt(:,1) + Urho(i,2)*phiValueInt(:,2) + Urho(i,3)*phiValueInt(:,3);
    Qrhou = Urhou(i,1)*phiValueInt(:,1) + Urhou(i,2)*phiValueInt(:,2) + Urhou(i,3)*phiValueInt(:,3);
    QE = UE(i,1)*phiValueInt(:,1) + UE(i,2)*phiValueInt(:,2) + UE(i,3)*phiValueInt(:,3);
    
    for i1 = 1:quad
        F = f(Qrho(i1),Qrhou(i1),QE(i1));
        f1 = F(1); f2 = F(2); f3 = F(3);
        for m = 1:k + 1
            dUrho(i,m) = dUrho(i,m) + hx1*weight(i1)*f1*phiValueIntx(i1,m);
            dUrhou(i,m) = dUrhou(i,m) + hx1*weight(i1)*f2*phiValueIntx(i1,m);
            dUE(i,m) = dUE(i,m) + hx1*weight(i1)*f3*phiValueIntx(i1,m);
        end
    end
    
    % 面积分 \hat -f(u_R-,u_R+)*v_R + \hat f(u_L-,u_L+)*v_L
    QrhoL = Urho(i,1)*phiValueL(1) + Urho(i,2)*phiValueL(2) + Urho(i,3)*phiValueL(3);
    QrhouL = Urhou(i,1)*phiValueL(1) + Urhou(i,2)*phiValueL(2) + Urhou(i,3)*phiValueL(3);
    QEL = UE(i,1)*phiValueL(1) + UE(i,2)*phiValueL(2) + UE(i,3)*phiValueL(3);
    
    QrhoLL = UrhoL(i,1)*phiValueR(1) + UrhoL(i,2)*phiValueR(2) + UrhoL(i,3)*phiValueR(3);
    QrhouLL = UrhouL(i,1)*phiValueR(1) + UrhouL(i,2)*phiValueR(2) + UrhouL(i,3)*phiValueR(3);
    QELL = UEL(i,1)*phiValueR(1) + UEL(i,2)*phiValueR(2) + UEL(i,3)*phiValueR(3);
    
    QrhoR = Urho(i,1)*phiValueR(1) + Urho(i,2)*phiValueR(2) + Urho(i,3)*phiValueR(3);
    QrhouR = Urhou(i,1)*phiValueR(1) + Urhou(i,2)*phiValueR(2) + Urhou(i,3)*phiValueR(3);
    QER = UE(i,1)*phiValueR(1) + UE(i,2)*phiValueR(2) + UE(i,3)*phiValueR(3);
    
    QrhoRR = UrhoR(i,1)*phiValueL(1) + UrhoR(i,2)*phiValueL(2) + UrhoR(i,3)*phiValueL(3);
    QrhouRR = UrhouR(i,1)*phiValueL(1) + UrhouR(i,2)*phiValueL(2) + UrhouR(i,3)*phiValueL(3);
    QERR = UER(i,1)*phiValueL(1) + UER(i,2)*phiValueL(2) + UER(i,3)*phiValueL(3);
    
    FR = f(QrhoR,QrhouR,QER); FRR = f(QrhoRR,QrhouRR,QERR);
    FL = f(QrhoL,QrhouL,QEL); FLL = f(QrhoLL,QrhouLL,QELL);
    
    [aR1,~,aR3] = eigenvalues(QrhoR,QrhouR,QER);
    [aRR1,~,aRR3] = eigenvalues(QrhoRR,QrhouRR,QERR);
    
    [aL1,~,aL3] = eigenvalues(QrhoL,QrhouL,QEL);
    [aLL1,~,aLL3] = eigenvalues(QrhoLL,QrhouLL,QELL);
    
    SRmax = max(aR1,aRR1); SRmin = min(aR3,aRR3);
    SLmax = max(aL1,aLL1); SLmin = min(aL3,aLL3);
    
    fhatR = flux(FR,FRR,QrhoR,QrhouR,QER,QrhoRR,QrhouRR,QERR,SRmax,SRmin);
    fhatL = flux(FLL,FL,QrhoLL,QrhouLL,QELL,QrhoL,QrhouL,QEL,SLmax,SLmin);
    
    fR1 = fhatR(1); fR2 = fhatR(2); fR3 = fhatR(3);
    fL1 = fhatL(1); fL2 = fhatL(2); fL3 = fhatL(3);
    
    for m = 1:k + 1
        dUrho(i,m) = dUrho(i,m) - fR1*phiValueR(m) + fL1*phiValueL(m);
        dUrhou(i,m) = dUrhou(i,m) - fR2*phiValueR(m) + fL2*phiValueL(m);
        dUE(i,m) = dUE(i,m) - fR3*phiValueR(m) + fL3*phiValueL(m);
    end
end

for m = 1:k + 1
    dUrho(:,m) = dUrho(:,m)/mm(m);
    dUrhou(:,m) = dUrhou(:,m)/mm(m);
    dUE(:,m) = dUE(:,m)/mm(m);
end

end

function y = f(rho,rhou,E)

global gamma

u = rhou/rho;

p = (gamma - 1)*(E - 0.5*rho*u^2);

y = [rho*u;rho*u^2 + p;u*(E + p)];

end
    
