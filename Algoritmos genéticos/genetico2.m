%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%ALGORITMOS GENETICOS, EJERCICIO DE DISTANCIAS
clear;

x = randi(255, 100, 1);
y = randi(255, 100, 1);
mCromosoma = zeros(100, 16);
nCromosoma = zeros(100, 16);
Pc = 0.75;  %0.25  mejor: 0.75
Pm = 0.01; %0.1  mejor: 0.01
desviacion = 10000;
minimo = 255000;

%% decimal a binario
cromosoma = [de2bi(x, 8) de2bi(y, 8)];
mitad = length(cromosoma(1,:))/2;
tam = mitad*2;
figure
hold on

%% binario a decimal escalado
for i=1:100 
   x(i) =  DecimalEscalado(cromosoma(i,1:mitad), 0, 12);
   y(i) =  DecimalEscalado(cromosoma(i,mitad+1:tam), 0, 12);
end

%% Funcion objetivo
f = 400*sqrt((x-5).^2 + (y-10).^2) + 300*sqrt((x-10).^2 + (y-5).^2) + 400*sqrt(x.^2 + (y-12).^2) + 600*sqrt((x-12).^2 + y.^2);


while (desviacion > 300) % minimo 300
    %% Fitness y proporcion del mismo
    fitness = 1./(1+f);
    tfit = sum(fitness);
    probmezcla = fitness./tfit;

    %Probabilidad acumulada
    probacum = cumsum(probmezcla);

    %% generar numeros aleatorios para conformar la nueva población
    %cruzada por ruleta
    randprob = rand(100, 1);
    mCromosoma = cromosoma;
    for i=1:100
        pre = find(probacum >= randprob(i));
        if ~isempty(pre)
            mCromosoma(i,:) = cromosoma(pre(1),:);
        end
    end

    %% hacer el cruce de cromosomas, el nuevo cromosoma reemplaza al padre #1
    probcruce = rand(100, 1);
    padres = find(probcruce < Pc);
    nCromosoma = mCromosoma;
    for i=1:length(padres)-1
        cant = randi(tam-1);
        nCromosoma(padres(i),:) = [mCromosoma(padres(i), 1:cant) mCromosoma(padres(i+1), cant+1:tam)];
    end
    nCromosoma(length(padres),:) = [mCromosoma(padres(length(padres)), 1:cant) mCromosoma(padres(1), cant+1:tam)];

    %% Calculo de la cantidad de mutaciones
    Tg = 100*tam;
    gM = Tg * Pm;

    %% Genes a mutar y su respectiva mutación
    genesM = randi(Tg, gM, 1);
    nCromosoma(genesM) = ~nCromosoma(genesM);
    cromosoma = nCromosoma;
    
    %% binario a decimal escalado
    for i=1:100 
       x(i) =  DecimalEscalado(cromosoma(i,1:mitad), 0, 12);
       y(i) =  DecimalEscalado(cromosoma(i,mitad+1:tam), 0, 12);
    end

    %% Funcion objetivo
    f = 400*sqrt((x-5).^2 + (y-10).^2) + 300*sqrt((x-10).^2 + (y-5).^2) + 400*sqrt(x.^2 + (y-12).^2) + 600*sqrt((x-12).^2 + y.^2);
    desviacion = std(f, 1);
end

%%  Resultados finales   
scatter3(x,y,f);
pos = find(min(f) == f);
scatter3(x(pos), y(pos), f(pos), 'r');
[x(pos) y(pos) f(pos)]
 