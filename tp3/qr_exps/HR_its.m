iteraciones=200;
load('HRs.mat');

h=figure;
plot([smooth(HR_100_0(1:iteraciones,1)),smooth(HR_100_0(1:iteraciones,2)),smooth(HR_100_0(1:iteraciones,3))]);
title('Hitrate a 500 imgs con 100 componentes, en funcion de iteraciones QR','FontSize',14.0);
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
legend('5k imagenes vecinos mas cercanos','15k imagenes vecinos mas cercanos','30k imagenes vecinos mas cercanos','Location','SouthEast');
axis([1 iteraciones 0.0 1]);
filename=sprintf('../informe/plots/HR_100_0');
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([smooth(HR_100_1(1:iteraciones,1)),smooth(HR_100_1(1:iteraciones,2)),smooth(HR_100_1(1:iteraciones,3))]);
axis([1 iteraciones 0.0 1])
legend('5k imagenes vecinos mas cercanos','15k imagenes vecinos mas cercanos','30k imagenes vecinos mas cercanos','Location','SouthEast');
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
title('Hitrate a 500 imgs con 100 componentes, en funcion de iteraciones QR','FontSize',14.0);
filename='../informe/plots/HR_100_1';
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([smooth(HR_10_0(1:iteraciones,1)),smooth(HR_10_0(1:iteraciones,2)),smooth(HR_10_0(1:iteraciones,3))]);
title('Hitrate a 500 imgs con 10 componentes, en funcion de iteraciones QR','FontSize',14.0);
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
legend('5k imagenes distancia a centro masas','15k imagenes distancia a centro masas','30k imagenes distancia a centro masas','Location','SouthEast');
axis([1 iteraciones 0.0 1]);
filename=sprintf('../informe/plots/HR_10_0');
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([smooth(HR_10_1(1:iteraciones,1)),smooth(HR_10_1(1:iteraciones,2)),smooth(HR_10_1(1:iteraciones,3))]);
title('Hitrate a 500 imgs con 10 componentes, en funcion de iteraciones QR','FontSize',14.0);
axis([1 iteraciones 0.0 1]);
legend('5k imagenes vecinos mas cercanos','15k imagenes vecinos mas cercanos','30k imagenes vecinos mas cercanos','Location','SouthEast');

xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
filename=sprintf('../informe/plots/HR_10_1');
saveas(h,filename,'png')
saveas(h,filename,'fig')


h=figure;
semilogy([SUM(:,1), SUM(:,2), SUM(:,3)]);
title('Suma elementos bajo la diagonal, en funcion de las iteraciones','FontSize',14.0);
legend('Suma 5k imagenes','Suma 15k imagenes','Suma 30k imagenes','Location','NorthEast');
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Suma valores bajo la diagonal','FontSize',14.0);
filename=sprintf('../informe/plots/SUM');
saveas(h,filename,'png')
saveas(h,filename,'fig')



h=figure;
semilogy([PROM(:,1), PROM(:,2), PROM(:,3)]);
title('Promedio elementos bajo la diagonal, en funcion de las iteraciones','FontSize',14.0);
legend('Promedio 5k imagenes','Promedio 15k imagenes','Promedio 30k imagenes','Location','NorthEast');
xlabel('Iteraciones de QR','FontSize',14.0);
ylabel('Promedio valores bajo la diagonal','FontSize',14.0);
filename=sprintf('../informe/plots/PROM');
saveas(h,filename,'png')
saveas(h,filename,'fig')
