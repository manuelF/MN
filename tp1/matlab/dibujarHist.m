function dibujarHist(data)

x_bins = min(data) : 0.02 : max(data);

% cuento cant repeticiones en rangos delimitados por x_bins
count = histc(data,x_bins);

bar(x_bins, count / sum(count));