% Plot USAT1
MTB = dlmread('MTB.42');
WBN = dlmread('WBN.42');
figure(1);clf;
subplot(2,1,1);
plot(MTB);hold on;grid on;title('MTB');
subplot(2,1,2);
plot(WBN);hold on;grid on;title('WBN');
