load('hr.mat')


h=figure;
plot([HR_5k(:,1), HR_15k(:,1), HR_30k(:,1)]);
legend('5k imagenes,a centros masa','15k imagenes, a centros masa','30k imagenes, a centros masa','Location','SouthEast');
axis([1 100 0 1]);
title('Hitrate usando deteccion a centros masa, en funcion cantidad componentes','FontSize',14.0);
xlabel('Cantidad de componentes principales tomadas','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
filename=sprintf('../informe/plots/HR_N2');
saveas(h,filename,'png')
saveas(h,filename,'fig')

h=figure;
plot([HR_5k(:,2), HR_15k(:,2), HR_30k(:,2)]);
legend('5k imagenes, vecinos cercanos','15k imagenes, vecinos cercanos','30k imagenes, vecinos cercanos','Location','SouthEast');
axis([1 100 0 1]);
title('Hitrate usando deteccion vecinos cercanos, en funcion cantidad componentes','FontSize',14.0);
xlabel('Cantidad de componentes principales tomadas','FontSize',14.0);
ylabel('Hitrate de deteccion','FontSize',14.0);
filename=sprintf('../informe/plots/HR_VEC');
saveas(h,filename,'png')
saveas(h,filename,'fig')