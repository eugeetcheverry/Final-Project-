% Plot TSAT
MTB = dlmread('MTB.42');
WBN = dlmread('wbn.42');
QBN = dlmread('qbn.42');
POSN = dlmread('PosN.42');
ILLUM = dlmread('Illum.42');
ALBEDO = dlmread('Albedo.42');
SVB = dlmread('svb.42');
figure(2);clf;
orbit=(1:length(POSN))/(96.77*60);xlabel('Orbits');
plot(orbit,POSN(:,3));hold on;grid on;title('POSN');
figure(1);clf;
subplot(3,1,1);
plot(orbit,MTB(:,1),'r');hold on;grid on;title('MTB');xlabel('Orbits');ylabel('Am^2');
plot(orbit,MTB(:,2),'g');hold on;grid on;title('MTB');xlabel('Orbits');ylabel('Am^2');
plot(orbit,MTB(:,3),'b');hold on;grid on;title('MTB');xlabel('Orbits');ylabel('Am^2');
axis([0 10 -0.03 0.03]);
subplot(3,1,2);
plot(orbit,WBN(:,1),'r');hold on;grid on;title('WBN');xlabel('Orbits');ylabel('rad/sec');
plot(orbit,WBN(:,2),'g');hold on;grid on;title('WBN');xlabel('Orbits');ylabel('rad/sec');
plot(orbit,WBN(:,3),'b');hold on;grid on;title('WBN');xlabel('Orbits');ylabel('rad/sec');
axis([0 10 -0.003 0.003]);
subplot(3,1,3);
plot(orbit,QBN(:,1),'r');hold on;grid on;title('QBN');
plot(orbit,QBN(:,2),'g');hold on;grid on;title('QBN');
plot(orbit,QBN(:,3),'b');hold on;grid on;title('QBN');
plot(orbit,QBN(:,4),'k');hold on;grid on;title('QBN');xlabel('Orbits');
axis([0 10 -1 1]);
figure(3);clf;
orbit=(1:length(POSN))/(96.77*60);xlabel('Orbits');
plot(orbit,QBN(:,1)*180./pi,'r');hold on;grid on;title('QBN');
plot(orbit,QBN(:,2)*180./pi,'g');hold on;grid on;title('QBN');
plot(orbit,QBN(:,3)*180./pi,'b');hold on;grid on;title('QBN');
title('QBN(0:2) [deg]');
figure(4);clf;
iiv=find(max(ILLUM')>0.);
plot(orbit(iiv),SVB(iiv,1),'r.');hold on;
plot(orbit(iiv),SVB(iiv,2),'g.');hold on;
plot(orbit(iiv),SVB(iiv,3),'b.');hold on;
grid on;title('SVB');
xlabel('Orbits');
figure(41);clf;
plot3(SVB(iiv,1),SVB(iiv,2),SVB(iiv,3),'k.');hold on;
figure(5);clf;
plot(orbit,ILLUM(:,1),'r');hold on;
plot(orbit,ILLUM(:,2),'r--');hold on;
plot(orbit,ILLUM(:,3),'g');hold on;
plot(orbit,ILLUM(:,4),'g--');hold on;
plot(orbit,ILLUM(:,5),'b');hold on;
plot(orbit,ILLUM(:,6),'b--');hold on;
%plot(orbit,ILLUM(:,7),'m');hold on;
%plot(orbit,ILLUM(:,8),'m--');hold on;
%plot(orbit,ILLUM(:,9),'c');hold on;
%plot(orbit,ILLUM(:,10),'c--');hold on;
%plot(orbit,ILLUM(:,11),'k');hold on;
%plot(orbit,ILLUM(:,12),'k--');hold on;

if 1==1
    
HWHL = dlmread('Hwhl.42');
figure(6);clf;
subplot(4,1,1);
plot(orbit,MTB(:,1),'r');hold on;grid on;title('MTB');xlabel('Orbits');ylabel('Am^2');
plot(orbit,MTB(:,2),'g');hold on;grid on;title('MTB');xlabel('Orbits');ylabel('Am^2');
plot(orbit,MTB(:,3),'b');hold on;grid on;title('MTB');xlabel('Orbits');ylabel('Am^2');
%axis([0 10 -0.01 0.01]);
subplot(4,1,2);
title('Wheel Speed');hold on;grid on;
plot(orbit,HWHL(:,1)/(2.1128e-6)/(2*pi)*60,'b');hold on;ylabel('RPM');grid on;
xlabel('Orbits');hold on;
%axis([0 10 -1500 1000]);
subplot(4,1,3);
plot(orbit,WBN(:,1),'r');hold on;grid on;title('WBN');xlabel('Orbits');ylabel('rad/sec');
plot(orbit,WBN(:,2),'g');hold on;grid on;title('WBN');xlabel('Orbits');ylabel('rad/sec');
plot(orbit,WBN(:,3),'b');hold on;grid on;title('WBN');xlabel('Orbits');ylabel('rad/sec');
subplot(4,1,4);
plot(orbit,QBN(:,1),'r');hold on;grid on;title('QBN');
plot(orbit,QBN(:,2),'g');hold on;grid on;title('QBN');
plot(orbit,QBN(:,3),'b');hold on;grid on;title('QBN');
plot(orbit,QBN(:,4),'k');hold on;grid on;title('QBN');xlabel('Orbits');
%axis([0 10 -1 1]);
%axis([0 10 -0.002 0.001]);
%legend('No Unloading','Unloading','Location','WestSouth')

end;

WF=6000;
ctita=zeros(size(SVB(:,1)));
ctita(iiv)=SVB(iiv,2);
cossunavg=conv(abs(ctita),ones(WF,1))/WF;
cossunavg=cossunavg(WF:(length(cossunavg)-WF));
figure(7);hold on;grid on;hold on;
title('Cos(TitaSol)');hold on;
plot((1:length(cossunavg))/(96.77*60),cossunavg,'b--');
xlabel('Orbits');
mean(cossunavg(10000:end))

RAAN = [90, 67.5, 45, 35, 30 , 22.5, 10, 0];
POWFACTOR = [0.57, 0.5845, 0.6418, 0.7373, 0.8018, 0.7954 0.8206 0.8248];
RAAN = [90, 67.5, 45, 35, 22.5, 10, 0];
POWFACTOR = [0.57, 0.5845, 0.6418, 0.7373, 0.7954 0.8206 0.8248];
RAAN =      [90,     67.5,   45,     35,  30,   22.5,  10,    0];
POWFACTOR = [0.6167, 0.6019, 0.6407, 0.7,   0.82614,0.8235,  0.8234,   0.8];
POWFACTORm =[0.46,   0.45,  0.474,    0.57,   0.65,     0.63,  0.65,     0.64];
RAANw =     [90,      67.5,   45,     35,     30,      22.5,    10,       0];
POWFACTORw =[0.6367, 0.661, 0.7452, 0.8395, 0.9938, 0.99864, 0.998936, 0.9983];
POWFACTORwm=[0.614 , 0.628, 0.71  , 0.7990, 0.96  , 0.996  , 0.9965,   0.985];
figure(10);clf;
title('Mean Power Generation Factor');hold on;
xlabel('RAAN');
grid on;
plot(RAAN,POWFACTOR,'ro-');hold on;
plot(RAAN,POWFACTORm,'r-.');hold on;
plot(RAANw,POWFACTORw,'bo-');hold on;
plot(RAANw,POWFACTORwm,'b-.');hold on;

p=polyfit(RAAN,POWFACTOR,4);
vRAAN=0:1:90;
vPOWFACTOR=polyval(p,vRAAN);
%plot(vRAAN,vPOWFACTOR,'r');hold on;
pw=polyfit(RAANw(1:(end-3)),POWFACTORw(1:(end-3)),2);
vRAANw=22:1:90;
vPOWFACTORw=polyval(pw,vRAANw);
%plot(vRAANw,vPOWFACTORw,'b');hold on;
plot(RAANw((end-2):end),POWFACTORw((end-2):end),'b-');hold on;
%legend('No Wheel','One Wheel','Fit No Wheel','Fit One Wheel','Location','NorthEast')


%chat;
%figure(100);clf;
%plot(orbit, vc, 'k');
%hold on;grid on;
%title('Estimation of the drag parameter');
%hold on;
%xlabel('Orbits');ylabel('Parameter');


%rref;
%figure(100);clf;
%plot(orbit, vr, 'k');
%hold on;grid on;
%title('Roll Error');
%hold on;
%xlabel('Orbits');ylabel('[degrees]');

if 1==0

mission;
figure(30);clf;
orbit2=(1:length(POSN))/(96.77*60);
plot(orbit2(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,1))*180/pi*2,'r');hold on;grid on;
plot(orbit2(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,2))*180/pi*2,'g');hold on;grid on;
plot(orbit2(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,3))*180/pi*2,'b');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,4))*180/pi*2,'k');hold on;grid on;
xlabel('[Seconds]');hold on;
ylabel('[Degrees]');hold on;
title('R Pointing Error [deg]');


figure(31);clf;
orbit2=(1:length(POSN))/(96.77*60);
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,1))*180/pi*2,'r');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,2))*180/pi*2,'g');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,3))*180/pi*2,'b');hold on;grid on;
plot(orbit2(length(POSN)-length(vm)+1:length(POSN)),(vm(:,5)),'k');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Nm]');hold on;
title('R Torque [Nm]');


%figure(30);clf;
%orbit2=(1:length(POSN))/(96.77*60);xlabel('Orbits');
%plot(orbit2(2:length(POSN)),(QBN(2:length(POSN),1)-QBN(1:(length(POSN)-1),1))*180/pi,'r');hold on;grid on;
%plot(orbit2(2:length(POSN)),(QBN(2:length(POSN),1)-QBN(1:(length(POSN)-1),2))*180/pi,'g');hold on;grid on;
%plot(orbit2(2:length(POSN)),(QBN(2:length(POSN),1)-QBN(1:(length(POSN)-1),3))*180/pi,'b');hold on;grid on;
%title('QBN(0:2) [deg]');

mission;
eMM = [];
for ii=1:120:length(POSN)
    MM = [[vm(ii,6) vm(ii,7) vm(ii,8)];[vm(ii,9) vm(ii,10) vm(ii,11)];[vm(ii,12) vm(ii,13) vm(ii,14)]];
    eMM = [eMM eig(MM)];
end;

figure(411);hold on;
plot((1:length(eMM))*2/0.962/100,eMM(1,:),'r');hold on;grid on;
plot((1:length(eMM))*2/0.962/100,eMM(2,:),'g');hold on;grid on;
plot((1:length(eMM))*2/0.962/100,eMM(3,:),'b');hold on;grid on;
title('Eigenvalues of \Gamma(t)');
xlabel('Orbits');

end;