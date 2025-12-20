% Plot TSAT
MTB = dlmread('MTB.42');
WBN = dlmread('wbn.42');
QBN = dlmread('qbn.42');
ILLUM = dlmread('Illum.42');
ALBEDO = dlmread('Albedo.42');
SVB = dlmread('svb.42');

VELN = dlmread('VelN.42');
POSN = dlmread('PosN.42');

orbit=(1:length(POSN))/(96.77*60);

if 1==1
vmean_elem=[];
voscul_elem=[];
figure(2);clf;
plot(orbit,POSN(:,3));hold on;grid on;title('POSN');
for io=1:length(orbit)
    [oscul_elem, mean_elem, ustin_elem] = testfastmean(POSN(io,:)',VELN(io,:)');
    vmean_elem=[vmean_elem mean_elem];
    voscul_elem=[voscul_elem oscul_elem];
end;
figure(22);hold on;
orbit=(1:length(vmean_elem))/(96.77*60);
subplot(3,2,1);plot(orbit,voscul_elem(1,:),'r');grid on;xlabel('Orbits');hold on;
subplot(3,2,1);plot(orbit,vmean_elem(1,:),'k');grid on;title('a [m]');xlabel('Orbits');
subplot(3,2,2);plot(orbit,voscul_elem(2,:),'r');grid on;title('e');xlabel('Orbits');hold on;
subplot(3,2,2);plot(orbit,vmean_elem(2,:),'k');grid on;title('e');
subplot(3,2,3);plot(orbit,voscul_elem(3,:),'r');grid on;title('i [deg]');xlabel('Orbits');hold on;
subplot(3,2,3);plot(orbit,vmean_elem(3,:)*180/pi,'k');grid on;title('i [deg]');
subplot(3,2,4);plot(orbit,voscul_elem(4,:),'r');grid on;title('RAAN [deg]');axis([0 max(orbit) 0 30]);xlabel('Orbits');hold on;
subplot(3,2,4);plot(orbit,vmean_elem(4,:)*180/pi,'k');grid on;title('RAAN [deg]');axis([0 max(orbit) 0 0.5]);
subplot(3,2,5);plot(orbit,vmean_elem(5,:)*180/pi,'k');grid on;title('w [deg]');xlabel('Orbits');
subplot(3,2,6);plot(orbit,vmean_elem(6,:)*180/pi,'k');grid on;title('nu [deg]');xlabel('Orbits');
end;

figure(1);clf;
subplot(3,1,1);
plot(orbit,MTB);hold on;grid on;title('MTB');xlabel('Orbits');hold on;ylabel('Am^2');grid on;hold on;
axis([0 2 -45 45]);
subplot(3,1,2);
plot(orbit,WBN(:,1)*180/pi,'r');hold on;grid on;title('WBN');xlabel('Orbits');
plot(orbit,WBN(:,2)*180/pi,'g');hold on;grid on;title('WBN');xlabel('Orbits');
plot(orbit,WBN(:,3)*180/pi,'b');hold on;grid on;title('WBN');xlabel('Orbits');
ylabel('deg/s');
axis([0 2 -2 2]);
%plot(orbit,sqrt(WBN(:,1).^2+WBN(:,2).^2+WBN(:,3).^2),'k');hold on;grid on;title('WBN');xlabel('Orbits');
subplot(3,1,3);
plot(orbit,QBN(:,1),'r');hold on;grid on;title('QBN');
plot(orbit,QBN(:,2),'g');hold on;grid on;title('QBN');
plot(orbit,QBN(:,3),'b');hold on;grid on;title('QBN');
plot(orbit,QBN(:,4),'k');hold on;grid on;title('QBN');xlabel('Orbits');
axis([0 2 -1 1]);
if 1==0
figure(3);clf;
orbit=(1:length(POSN))/(96.77*60);xlabel('Orbits');
plot(orbit,QBN(:,1)*180./pi,'r');hold on;grid on;title('QBN');
plot(orbit,QBN(:,2)*180./pi,'g');hold on;grid on;title('QBN');
plot(orbit,QBN(:,3)*180./pi,'b');hold on;grid on;title('QBN');
title('QBN(0:2) [deg]');
end;
figure(4);clf;
iiv=find(max(ILLUM')>0.);
plot(orbit(iiv),SVB(iiv,1),'r.');hold on;
plot(orbit(iiv),SVB(iiv,2),'g.');hold on;
plot(orbit(iiv),SVB(iiv,3),'b.');hold on;
grid on;title('SVB');
xlabel('Orbits');
axis([0 2 -1 1]);
if 1==0
figure(41);clf;
plot3(SVB(iiv,1),SVB(iiv,2),SVB(iiv,3),'k.');hold on;
end;
figure(5);clf;
plot(orbit,ILLUM(:,1),'r');hold on;
plot(orbit,ILLUM(:,2),'r--');hold on;
plot(orbit,ILLUM(:,3),'g');hold on;
plot(orbit,ILLUM(:,4),'g--');hold on;
plot(orbit,ILLUM(:,5),'b');hold on;
plot(orbit,ILLUM(:,6),'b--');hold on;
plot(orbit,ILLUM(:,7),'m');hold on;
plot(orbit,ILLUM(:,8),'m--');hold on;
plot(orbit,ILLUM(:,9),'c');hold on;
plot(orbit,ILLUM(:,10),'c--');hold on;
%plot(orbit,ILLUM(:,11),'k');hold on;
%plot(orbit,ILLUM(:,12),'k--');hold on;
grid on;title('CSS Illumination');
xlabel('Orbits');
HWHL = dlmread('Hwhl.42');
figure(6);hold on;
orbit=(1:length(POSN))/(96.77*60);xlabel('Orbits');
title('Wheel Momentum Nms');hold on;grid on;
%plot(orbit,HWHL(:,1)/(3.5e-3)/(2*pi)*60,'r--');hold on;
%plot(orbit,HWHL(:,2)/(3.5e-3)/(2*pi)*60,'g--');hold on;
%plot(orbit,HWHL(:,3)/(3.5e-3)/(2*pi)*60,'b--');hold on;
plot(orbit,HWHL(:,1),'r');hold on;
plot(orbit,HWHL(:,2),'g');hold on;
plot(orbit,HWHL(:,3),'b');hold on;
xlabel('Orbits');hold on;
axis([0 2 -0.6 0.3]);
ylabel('[Nms]');
%legend('No Unloading','Unloading','Location','WestSouth')

if 1==0

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

end;

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


missionf;
orbit=0.5*(1:length(vm))/(96.77*60);

if 1==1

    figure(33);clf;
    subplot(3,2,1);plot(orbit,vm(:,28));grid on;title('a');
    subplot(3,2,2);plot(orbit,vm(:,29));grid on;title('e');
    subplot(3,2,3);plot(orbit,vm(:,30)*180/pi);grid on;title('i [deg]');
    subplot(3,2,4);plot(orbit,vm(:,31)*180/pi);grid on;title('RAAN [deg]');
    subplot(3,2,5);plot(orbit,vm(:,32)*180/pi);grid on;title('w [deg]');
    subplot(3,2,6);plot(orbit,vm(:,34)*180/pi);grid on;title('lm [deg]');

    figure(34);clf;
    subplot(3,1,1);
    title('RTN Thrust [N]');hold on;grid on;
    plot(orbit,vm(:,35),'r');hold on;
    plot(orbit,vm(:,36),'g');hold on;
    plot(orbit,vm(:,37),'b');hold on;
    xlabel('Orbits');
    subplot(3,1,2);
    title('DeltaV [s]');hold on;grid on;
    timp=0.5*cumsum(sqrt(vm(:,35).^2+vm(:,36).^2+vm(:,37).^2))/133;
    plot(orbit,timp,'k');hold on;
    dv400=0.11*ones(size(vm(:,1)));
    plot(orbit,dv400,'k--');hold on,
    xlabel('Orbits');
    subplot(3,1,3);
    plot(orbit,vm(:,28));grid on;title('a [m]');
    xlabel('Orbits');

    figure(35);
    title('Total Impulse [Ns]');hold on;grid on;
    timp=0.5*cumsum(sqrt(vm(:,35).^2+vm(:,36).^2+vm(:,37).^2));
    plot(orbit,timp,'k');hold on;
    xlabel('Orbits');

    figure(36);
    title('DeltaV [s]');hold on;grid on;
    timp=0.5*cumsum(sqrt(vm(:,35).^2+vm(:,36).^2+vm(:,37).^2))/133;
    plot(orbit,timp,'k');hold on;
    xlabel('Orbits');

end;

if 1==0
   
figure(30);clf;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,1))*180/pi*2,'r');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,2))*180/pi*2,'g');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,3))*180/pi*2,'b');hold on;grid on;
r1=8000:8500;
subplot(1,3,1);
plot(orbit(r1),asin(vm(r1     ,1))*180/pi*2,'k');hold on;grid on;
%plot((vm(r1+2500,1)),'k');hold on;grid on;
%plot((vm(r1+2500,1)),'k');hold on;grid on;
%plot((vm(r1+5000,1)),'k');hold on;grid on;
%plot((vm(r1+7500,1)),'k');hold on;grid on;
%plot((vm(r1+10000,1)),'k');hold on;grid on;
%plot((vm(r1+12500,1)),'k');hold on;grid on;
xlabel('[Seconds]');hold on;
ylabel('[Degrees]');hold on;
title('T Pointing Error [deg]');
subplot(1,3,2);
plot(orbit(r1),asin(vm(r1     ,2))*180/pi*2,'k');hold on;grid on;
%plot((vm(r1+2500,2)),'k');hold on;grid on;
%plot((vm(r1+2500,2)),'k');hold on;grid on;
%plot((vm(r1+5000,2)),'k');hold on;grid on;
%plot((vm(r1+7500,2)),'k');hold on;grid on;
%plot((vm(r1+10000,2)),'k');hold on;grid on;
%plot((vm(r1+12500,2)),'k');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Degrees]');hold on;
title('N Pointing Error [deg]');
subplot(1,3,3);
plot(orbit(r1),asin(vm(r1     ,3))*180/pi*2,'k');hold on;grid on;
%plot((vm(r1+2500,3)),'k');hold on;grid on;
%plot((vm(r1+2500,3)),'k');hold on;grid on;
%plot((vm(r1+5000,3)),'k');hold on;grid on;
%plot((vm(r1+7500,3)),'k');hold on;grid on;
%plot((vm(r1+10000,3)),'k');hold on;grid on;
%plot((vm(r1+12500,3)),'k');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,4))*180/pi*2,'k');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Degrees]');hold on;
title('R Pointing Error [deg]');

end;

if 1==0
figure(29);clf;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,1))*180/pi*2,'r');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,2))*180/pi*2,'g');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,3))*180/pi*2,'b');hold on;grid on;
r1=1:length(vm);
plot(orbit(r1),asin(vm(r1     ,1))*180/pi*2,'r');hold on;grid on;
plot(orbit(r1),asin(vm(r1     ,2))*180/pi*2,'g');hold on;grid on;
plot(orbit(r1),asin(vm(r1     ,3))*180/pi*2,'b');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Degrees]');hold on;
title('TNR Pointing Error [deg]');
end;


figure(292);clf;
r1=1:length(vm);
plot(r1/2,vm(r1,13),'r');hold on;grid on;
plot(r1/2,vm(r1,14),'g');hold on;grid on;
plot(r1/2,vm(r1,15),'b');hold on;grid on;
plot(r1/2,vm(r1,16),'k');hold on;grid on;
xlabel('[Seconds]');hold on;
ylabel('[Degrees]');hold on;
title('QLN [deg]');



figure(293);clf;
r1=1:length(vm);
plot(r1/2,vm(r1,17),'r');hold on;grid on;
plot(r1/2,vm(r1,18),'g');hold on;grid on;
plot(r1/2,vm(r1,19),'b');hold on;grid on;
plot(r1/2,vm(r1,20),'k');hold on;grid on;
xlabel('[Seconds]');hold on;
ylabel('[Degrees]');hold on;
title('qbn');


figure(294);clf;
r1=1:length(vm);
title('chat x[r], chaty[g], chatz[b]');
plot(r1/2,vm(r1,21),'r');hold on;grid on;
plot(r1/2,vm(r1,22),'g');hold on;grid on;
plot(r1/2,vm(r1,23),'b');hold on;grid on;
xlabel('[Seconds]');hold on;
ylabel('[chat]');hold on;


figure(291);clf;
r1=1:length(vm);
plot(r1/2-6000*2,asin(vm(r1     ,1))*180/pi*2,'r');hold on;grid on;
plot(r1/2-6000*2,asin(vm(r1     ,2))*180/pi*2,'g');hold on;grid on;
plot(r1/2-6000*2,asin(vm(r1     ,3))*180/pi*2,'b');hold on;grid on;
%plot(r1/2,vm(r1     ,24)*180/pi+60,'m');hold on;grid on;
%plot(r1/2,vm(r1     ,25)*180/pi+60,'c');hold on;grid on;
xlabel('[Seconds]');hold on;
ylabel('[Degrees]');hold on;
axis([0 600 -0.1 0.1]);
title('SAR Pointing Error [deg]');


figure(2910);clf;
plot(r1/2-6000*2-600,asin(vm(r1     ,1))*180/pi*2+vm(r1     ,25)*180/pi,'r');hold on;
plot(r1/2-6000*2-600,-vm(r1     ,24)*180/pi,'b');hold on;
plot(r1/2-6000*2-600, vm(r1     ,25)*180/pi,'g');hold on;
plot(r1/2-6000*2-600,-vm(r1     ,27)*180/pi,'k');hold on;
grid on;
xlabel('[Seconds]');hold on;
ylabel('[Degrees]');hold on;
title('SAR Pointing (r), filter 1 (b), filter 2 (b) and reference (k) [deg]');

figure(31);clf;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,1))*180/pi*2,'r');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,2))*180/pi*2,'g');hold on;grid on;
%plot(orbit(length(POSN)-length(vm)+1:length(POSN)),asin(vm(:,3))*180/pi*2,'b');hold on;grid on;
plot(r1/2-6000,(vm(:,5)),'r');hold on;grid on;
plot(r1/2-6000,(vm(:,6)),'g');hold on;grid on;
plot(r1/2-6000,(vm(:,7)),'b');hold on;grid on;
xlabel('[Seconds]');hold on;
axis([0 600 -0.1 0.1]);
ylabel('[Nm]');hold on;
title('RW Torque [Nm]');


if 1==0

figure(32);clf;
plot(orbit,(vm(:,8))*180/pi,'k');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Degrees]');hold on;
title('Phi [Nm]');


figure(33);clf;
plot(orbit,(vm(:,9)),'r');hold on;grid on;
plot(orbit,(vm(:,10)),'k');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Seconds]');hold on;
title('secondsmis [r] secondsnmis [k]');


%figure(30);clf;
%orbit=(1:length(POSN))/(96.77*60);xlabel('Orbits');
%plot(orbit(2:length(POSN)),(QBN(2:length(POSN),1)-QBN(1:(length(POSN)-1),1))*180/pi,'r');hold on;grid on;
%plot(orbit(2:length(POSN)),(QBN(2:length(POSN),1)-QBN(1:(length(POSN)-1),2))*180/pi,'g');hold on;grid on;
%plot(orbit(2:length(POSN)),(QBN(2:length(POSN),1)-QBN(1:(length(POSN)-1),3))*180/pi,'b');hold on;grid on;
%title('QBN(0:2) [deg]');

figure(33);clf;
plot(orbit,(vm(:,11))*180/pi,'r');hold on;grid on;
plot(orbit,(vm(:,12))*180/pi,'k');hold on;grid on;
plot(orbit,cumsum((vm(:,12)))*180/pi/2,'k--');hold on;grid on;
%plot(orbit(2:end),(vm(2:end,11))*180/pi-(vm(1:(end-1),11))*180/pi,'c--');hold on;grid on;
xlabel('[Orbits]');hold on;
ylabel('[Degrees]');hold on;
title('Phi [r] dPhi [k]');

end;

if 1==0
    
% Floquet
veps = [0.0001 0.001 0.01; 0.1; 1];

%X0=[];XT=[];
WBN = dlmread('wbn.42');
WBNS3 = WBN;
WBN = WBNS3;
NT = min(find(orbit>=1));
XT = [XT WBN(NT,:)'];
X0 = [X0 WBN(1,:)'];

B = inv(X0)*XT;
mul_floquet = eig(B);
exp_floquet = log(mul_floquet)/(NT/2);

vexp_floquet = [vexp_floquet exp_floquet];

save("vexp_floquet.mat",vexp_floquet,veps);

figure(1000);clf;
title('Maximum Floquet Exponents (real part) for k_v=10 and SSO');hold on;
plot(veps,max(real(vexp_floquet)),'r');hold on;
xlabel('\epsilon');
grid on;



end;