
%PLOT_CURRENTS Summary of this function goes here

t = time; 
ax1 = subplot(4,4,1),hold on, plot(t,y(:,39),'LineWidth',2); ylabel('Em (mV)');

ax2 = subplot(4,4,2),hold on,plot(tArray, Ica(:),'LineWidth',2); ylabel('ICAL (pA/pF)');

ax3 = subplot(4,4,3),hold on,plot(t,100*y(:,87+45+2)./LCCtotDyad,'LineWidth',2); ylabel('LTCCp (%)');

ax4 = subplot(4,4,4),hold on,plot(t,y(:,38),'LineWidth',2); ylabel('[Ca]i (mM)');

ax5 = subplot(4,4,5),hold on,plot(t,y(:,30)+y(:,31),'LineWidth',2); ylabel('tot [Ca]SR (mM)'); xlabel('Time (s)');

ax6 = subplot(4,4,6),hold on, plot(t, y(:,31),'LineWidth',2); ylabel('Free SR Ca(mM)'); xlabel('Time (s)');

xx = find(tArray>60,1,'first'); xxx=find(tArray>1000,1,'first'); 
ax7 = subplot(4,4,7),hold on, plot(tArray(1:xx), Jleak(1:xx,1),'LineWidth',2); ylabel('RyR Jleak Systolic (pA/pF)'); xlabel('Time (s)');

ax8 = subplot(4,4,8),hold on, plot(tArray(xxx:end), Jleak(xxx:end,1),'LineWidth',2); ylabel('RyR Jleak Diastolic (pA/pF)'); xlabel('Time (s)');


ax9 = subplot(4,4,9),hold on, plot(tArray, Jserca(:,:),'LineWidth',2); ylabel('JSerca'); xlabel('Time (s)');

ax10 = subplot(4,4,10),hold on, plot(t, y(:,46),'LineWidth',2); ylabel('Background Ca(mM)'); xlabel('Time (s)');

ax11 = subplot(4,4,11),hold on, plot(tArray, Incx(:,:),'LineWidth',2); ylabel('Incx'); xlabel('Time (s)');

ax12 = subplot(4,4,12),hold on, plot(t, y(:,36),'LineWidth',2); ylabel('Cleft Ca(mM)'); xlabel('Time (s)');

ax13 = subplot(4,4,13), plot(t,y(:,14),'LineWidth',2); ylabel('RyR R'), hold on

ax14 = subplot(4,4,14), plot(t,y(:,15),'LineWidth',2); ylabel('RyR O'), hold on

ax15 = subplot(4,4,15), plot(t,y(:,16),'LineWidth',2); ylabel('RyR I'), hold on

ax16 = subplot(4,4,16), plot(tArray,Jrel,'LineWidth',2); ylabel('Jrel'), hold on

linkaxes([ax1  ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16],'x')


