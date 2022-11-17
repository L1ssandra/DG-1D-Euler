% flash1Dsystem.m
h = figure();
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

Uflashp = (gamma - 1)*(UflashE - 0.5*Uflashrhou.^2./Uflashrho);
Uflashu = Uflashrhou./Uflashrho;

Uplot = Uflashp;

xlabel('X');
ylabel('Q');
d = max(max(Uplot)) - min(min(Uplot));
s = [Xc(1),Xc(end),min(min(Uplot)) - 0.05*d,max(max(Uplot)) + 0.05*d];
%s = [Xc(1),Xc(end),0,10];
TT = Frame; % 帧数
t0 = T(end)/TT;

P = 1; % 目标分量

for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
   plot(Xc,Uplot(:,j),'b-');
   axis(s);
   pause(0.001);
end

% plot(Xc,QF(:,1,j),'r*',Xc,QR(:,1),'b-');
% axis(s);
% pause(0.001);