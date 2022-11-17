% test.m 测试三种Flux的精度

global type

type = 3;
Nx = 2500;
main
XcReal = Xc;
UflashrhoReal = Uflashrho;
UflashrhouReal = Uflashrhou;
Uflashp = (gamma - 1)*(UflashE - 0.5*Uflashrhou.^2./Uflashrho);
UflashuReal = UflashrhouReal./UflashrhoReal;
UflashpReal = Uflashp;
type = 1;
Nx = 200;
main
Uflashrho1 = Uflashrho;
Uflashrhou1 = Uflashrhou;
Uflashp = (gamma - 1)*(UflashE - 0.5*Uflashrhou.^2./Uflashrho);
Uflashu1 = Uflashrhou1./Uflashrho1;
Uflashp1 = Uflashp;
type = 2;
Nx = 200;
main
Uflashrho2 = Uflashrho;
Uflashrhou2 = Uflashrhou;
Uflashp = (gamma - 1)*(UflashE - 0.5*Uflashrhou.^2./Uflashrho);
Uflashu2 = Uflashrhou2./Uflashrho2;
Uflashp2 = Uflashp;
type = 3;
Nx = 200;
main
Uflashrho3 = Uflashrho;
Uflashrhou3 = Uflashrhou;
Uflashp = (gamma - 1)*(UflashE - 0.5*Uflashrhou.^2./Uflashrho);
Uflashu3 = Uflashrhou3./Uflashrho3;
Uflashp3 = Uflashp;

figure(1)
d = max(max(UflashrhoReal)) - min(min(UflashrhoReal));
s = [Xc(1),Xc(end),min(min(UflashrhoReal)) - 0.05*d,max(max(UflashrhoReal)) + 0.05*d];
ylabel('density')
plot(Xc,Uflashrho1(:,end),'r.',Xc,Uflashrho2(:,end),'m.',Xc,Uflashrho3(:,end),'c.',XcReal,UflashrhoReal(:,end),'b-');
legend('L-F','HLL','HLLC','exact')
axis(s)

figure(2)
d = max(max(UflashuReal)) - min(min(UflashuReal));
s = [Xc(1),Xc(end),min(min(UflashuReal)) - 0.05*d,max(max(UflashuReal)) + 0.05*d];
ylabel('Velocity')
plot(Xc,Uflashu1(:,end),'r.',Xc,Uflashu2(:,end),'m.',Xc,Uflashu3(:,end),'c.',XcReal,UflashuReal(:,end),'b-');
legend('L-F','HLL','HLLC','exact')
axis(s)

figure(3)
d = max(max(UflashpReal)) - min(min(UflashpReal));
s = [Xc(1),Xc(end),min(min(UflashpReal)) - 0.05*d,max(max(UflashpReal)) + 0.05*d];
ylabel('Pressure')
plot(Xc,Uflashp1(:,end),'r.',Xc,Uflashp2(:,end),'m.',Xc,Uflashp3(:,end),'c.',XcReal,UflashpReal(:,end),'b-');
legend('L-F','HLL','HLLC','exact')
axis(s)

