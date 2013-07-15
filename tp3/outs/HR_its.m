iteraciones=200
h=figure;
plot([smooth(HR_100_0(1:iteraciones,1)),smooth(HR_100_0(1:iteraciones,2)),smooth(HR_100_0(1:iteraciones,3))]);
title('Hitrate a 500 imgs con 100 componentes, en funcion de iteraciones QR','FontSize',14.0);
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
legend('5k VecinosCercanos','15k VecinosCercanos','30k VecinosCercanos','Location','SouthEast');
axis([1 iteraciones 0.0 1]);
filename=sprintf('../informe/plots/HR_100_0');
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([smooth(HR_100_1(1:iteraciones,1)),smooth(HR_100_1(1:iteraciones,2)),smooth(HR_100_1(1:iteraciones,3))]);
axis([1 iteraciones 0.0 1])
legend('5k DistPromedios','15k DistPromedios','30k DistPromedios','Location','SouthEast');
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
title('Hitrate a 500 imgs con 100 componentes, en funcion de iteraciones QR','FontSize',14.0);
filename='../informe/plots/HR_100_1';
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([smooth(HR_10_1(1:iteraciones,1)),smooth(HR_10_1(1:iteraciones,2)),smooth(HR_10_1(1:iteraciones,3))]);
title('Hitrate a 500 imgs con 10 componentes, en funcion de iteraciones QR','FontSize',14.0);
axis([1 iteraciones 0.0 1]);
legend('5k DistPromedios','15k DistPromedios','30k DistPromedios','Location','SouthEast');
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
filename=sprintf('../informe/plots/HR_10_1');
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([smooth(HR_10_0(1:iteraciones,1)),smooth(HR_10_0(1:iteraciones,2)),smooth(HR_10_0(1:iteraciones,3))]);
title('Hitrate a 500 imgs con 10 componentes, en funcion de iteraciones QR','FontSize',14.0);
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
legend('5k VecinosCercanos','15k VecinosCercanos','30k VecinosCercanos','Location','SouthEast');
axis([1 iteraciones 0.0 1]);
filename=sprintf('../informe/plots/HR_10_0');
saveas(h,filename,'png')
saveas(h,filename,'fig')