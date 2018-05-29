%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%ALGORITMOS GENETICOS

%% Esta función convierte un numero de binario a decimal dentro de la escala dada
function decimalescalado = DecimalEscalado(binarios, lower, upper)

k = 1:length(binarios);
total = sum(2.^(k-1));
num = binarios*(2.^(k-1))';
decimalescalado = lower + num*(upper - lower)/ total;