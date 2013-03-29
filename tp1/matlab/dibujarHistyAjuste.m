function dibujarHistyAjuste(data, sigma,beta,lambda)
%
% Dado un vector de datos data y los parametros sigma, beta y lambda de una
% distribucion Gamma Generalizada, se grafica el histograma de los datos y
% la funcion de densidad de la DGG.
%


x_bins = min(data) : 0.02 : max(data);

% cuento cant repeticiones en rangos delimitados por x_bins
count = histc(data,x_bins);

figure;
bar(x_bins, count / sum(count));

hold on
y = GGDpdf_c(x_bins,sigma,beta,lambda);
plot(x_bins, y/sum(y),'r')

