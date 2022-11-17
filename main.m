% main.m
%clear;clc

global hx hx1 Nx k dim gamma quad bc
xa = 0;
xb = 1;
tend = 0.012;
k = 2;
RKorder = 1;
dim = 3;
%Nx = 60;
gamma = 1.4;
quad = 2;

global type
%type = 3;

bc = 2;

if RKorder == 1
    CFL = 0.05;
elseif RKorder == 3
    CFL = 0.1;
end


hx = (xb - xa)/Nx;
hx1 = hx/2;

X = xa:hx:xb;
Xc = (X(1:end - 1) + X(2:end))/2;

% rho = @(x) 1 + 0.2*sin(pi*x);
% u = @(x) 1 + 0.*x;
% p = @(x) 1 + 0.*x;

rho = @(x) rho0(x);
u = @(x) u0(x);
p = @(x) p0(x);

U1 = @(x) rho(x);
U2 = @(x) rho(x).*u(x);
U3 = @(x) p(x)./(gamma - 1) + 0.5.*rho(x).*u(x).*u(x);

Urho = zeros(Nx,k + 1);
Urhou = zeros(Nx,k + 1);
UE = zeros(Nx,k + 1);

global phiValueInt phiValueIntx phiValueL phiValueR weight mm

% 基函数
phi1 = @(x) 1 + 0.*x;
phi2 = @(x) x./hx1;
phi3 = @(x) (x./hx1).^2 - 1/3;

phi1x = @(x) 0.*x;
phi2x = @(x) 1/hx1 + 0.*x;
phi3x = @(x) 2.*x./(hx1^2);

% Gauss积分点
if quad == 2
    lambda(1) = 0.5773502691896257645091488;
    lambda(2) = -0.5773502691896257645091488;
    weight(1) = 1;
    weight(2) = 1;
elseif quad == 4
    lambda(1) = -0.8611363115940525752239465;     weight(1) = 0.3478548451374538573730639;
    lambda(2) = -0.3399810435848562648026658;     weight(2) = 0.6521451548625461426269361;
    lambda(3) = 0.3399810435848562648026658;     weight(3) = 0.6521451548625461426269361;
    lambda(4) = 0.8611363115940525752239465;     weight(4) = 0.3478548451374538573730639;
end

lambdax = hx1*lambda;

% 储存基函数的值
phiValueInt = zeros(quad,k + 1);
for i = 1:quad
    phiValueInt(i,1) = phi1(lambdax(i));
    phiValueInt(i,2) = phi2(lambdax(i));
    phiValueInt(i,3) = phi3(lambdax(i));
    
    phiValueIntx(i,1) = phi1x(lambdax(i));
    phiValueIntx(i,2) = phi2x(lambdax(i));
    phiValueIntx(i,3) = phi3x(lambdax(i));
end
phiValueR = zeros(1,k + 1);
phiValueL = zeros(1,k + 1);

phiValueR(1) = phi1(hx1);
phiValueR(2) = phi2(hx1);
phiValueR(3) = phi3(hx1);

phiValueL(1) = phi1(-hx1);
phiValueL(2) = phi2(-hx1);
phiValueL(3) = phi3(-hx1);

global phiValue0

phiValue0 = zeros(1,k + 1);
phiValue0(1) = phi1(0);
phiValue0(2) = phi2(0);
phiValue0(3) = phi3(0);

% 质量
mm = [hx,hx/3,4*hx/45];

% L2投影
for i = 1:Nx
    for m = 1:k + 1
        for i1 = 1:quad
            
            Urho(i,m) = Urho(i,m) + hx1*weight(i1)*U1(Xc(i) + lambdax(i1))*phiValueInt(i1,m)/mm(m);
            Urhou(i,m) = Urhou(i,m) + hx1*weight(i1)*U2(Xc(i) + lambdax(i1))*phiValueInt(i1,m)/mm(m);
            UE(i,m) = UE(i,m) + hx1*weight(i1)*U3(Xc(i) + lambdax(i1))*phiValueInt(i1,m)/mm(m);
            
        end
    end
end

t = 0; T = 0;
Uflashrho = Urho(:,1);
Uflashrhou = Urhou(:,1);
UflashE = UE(:,1);

Frame = 400;
t1 = tend/Frame;
count = 1;

while t < tend
    
    alphamax = 1;
    for i = 1:Nx
        [alpha1,~,alpha2] = eigenvalues(Urho(i,1),Urhou(i,1),UE(i,1));
        if max(abs(alpha1),abs(alpha2)) > alphamax
            alphamax = max(abs(alpha1),abs(alpha2));
        end
    end
    
    dt = CFL/(alphamax/hx);
    
    if t + dt < tend
        t = t + dt;
    else
        dt = tend - t;
        t = tend;
    end
    
    
    if RKorder == 1
        [dUrho,dUrhou,dUE] = Lh(Urho,Urhou,UE);
        Urho = Urho + dt*dUrho;
        Urhou = Urhou + dt*dUrhou;
        UE = UE + dt*dUE;
        if k > 0
            [Urho,Urhou,UE] = TVB_Limiter(Urho,Urhou,UE);
            [Urho,Urhou,UE] = pp_Limiter(Urho,Urhou,UE);
        end
    elseif RKorder == 3
        % Stage 1
        [dUrho,dUrhou,dUE] = Lh(Urho,Urhou,UE);
        Urho1 = Urho + dt*dUrho;
        Urhou1 = Urhou + dt*dUrhou;
        UE1 = UE + dt*dUE;
        if k > 0
            [Urho1,Urhou1,UE1] = TVB_Limiter(Urho1,Urhou1,UE1);
            [Urho1,Urhou1,UE1] = pp_Limiter(Urho1,Urhou1,UE1);
        end
        % Stage 2
        [dUrho,dUrhou,dUE] = Lh(Urho1,Urhou1,UE1);
        Urho2 = (3/4)*Urho + (1/4)*Urho1 + (1/4)*dt*dUrho;
        Urhou2 = (3/4)*Urhou + (1/4)*Urhou1 + (1/4)*dt*dUrhou;
        UE2 = (3/4)*UE + (1/4)*UE1 + (1/4)*dt*dUE;
        if k > 0
            [Urho2,Urhou2,UE2] = TVB_Limiter(Urho2,Urhou2,UE2);
            [Urho2,Urhou2,UE2] = pp_Limiter(Urho2,Urhou2,UE2);
        end
        % Stage 3
        [dUrho,dUrhou,dUE] = Lh(Urho2,Urhou2,UE2);
        Urho = (1/3)*Urho + (2/3)*Urho2 + (2/3)*dt*dUrho;
        Urhou = (1/3)*Urhou + (2/3)*Urhou2 + (2/3)*dt*dUrhou;
        UE = (1/3)*UE + (2/3)*UE2 + (2/3)*dt*dUE;
        if k > 0
            [Urho,Urhou,UE] = TVB_Limiter(Urho,Urhou,UE);
            [Urho,Urhou,UE] = pp_Limiter(Urho,Urhou,UE);
        end
    end
    
    if t >= count*t1
        T = [T;t];
        Uflashrho(:,end + 1) = Urho(:,1);
        Uflashrhou(:,end + 1) = Urhou(:,1);
        UflashE(:,end + 1) = UE(:,1);
    end
    
    if dt < 1e-9
        t = tend;
    end
    
    fprintf('已经迭代到 t = %d\n',t)
    
end

Equad = zeros(quad*Nx,1);
for i = 1:Nx
    Equad((i - 1)*quad + 1:i*quad) = Urho(i,1)*phiValueInt(:,1) + Urho(i,2)*phiValueInt(:,2) + Urho(i,3)*phiValueInt(:,3) - U1((Xc(i) + lambdax)');
end

L8 = max(abs(Equad));

%flash1Dsystem