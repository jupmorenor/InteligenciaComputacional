%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%ALGORITMOS GENETICOS, EJERCICIO DE MULTIPLES OBJETIVOS 
%CRITERIO DE PARETO
clear;
tic
x = randi(255, 100, 1);
y = randi(255, 100, 1);
z = randi(255, 100, 1);
x1 = zeros(100, 1);
y1 = zeros(100, 1);
z1 = zeros(100, 1);
x2 = zeros(100, 1);
y2 = zeros(100, 1);
z2 = zeros(100, 1);
mCromosoma1 = zeros(100, 24);
mCromosoma2 = zeros(100, 24);
nCromosoma1 = zeros(100, 24);
nCromosoma2 = zeros(100, 24);
fop = zeros(0, 3);
Pc = 0.25;
Pm = 0.1;
r = [1; 1];
desviacion1 = 10000;
desviacion2 = 10000;
counter = 1;

%% Cromosoma inicial decimal a binario
cromosoma = [de2bi(x, 8) de2bi(y, 8) de2bi(z, 8)];
mitad = length(cromosoma(1,:))/3;
tam = mitad*3;

%% binario a decimal escalado
for c=1:100 
   x(c) = DecimalEscalado(cromosoma(c,1:mitad), 0, 32.5);
   y(c) = DecimalEscalado(cromosoma(c,mitad+1:mitad*2), 0, 43.3);
   z(c) = DecimalEscalado(cromosoma(c,mitad*2+1:tam), 0, 50);
end

%% PROBLEMA A SOLUCIONAR
%Funcion objetivo 1 a maximizar
f1 = 10*x + 9*y + 8*z;
%Funcion objetivo 2 a minimizar
f2 = 10*x + 6*y + 3*z;
%% Restricciones globales 
g1 = 130 - 4*x - 3*y - 2*z;
g2 = 100 - 3*x - 2*y - 2*z;

%while (desviacion1 > 100 && desviacion2 > 75) 
for l=1:5000
    %% Satisfaccion de las restricciones
    neg1 = find(g1 >= 0);
    neg2 = find(g2 >= 0);
    
    g1(neg1) = 0;
    g2(neg2) = 0; 
    
    %% Frente optimo de pareto
    fop1 = FrenteOptimo(f1, f2, g1, g2);
    if fop1(1) > 0
        fop = cat(1, fop, [x(fop1) y(fop1) z(fop1)]);
    end
    
    fpareto1 = 10*fop(:, 1) + 9*fop(:, 2) + 8*fop(:, 3);
    fpareto2 = 10*fop(:, 1) + 6*fop(:, 2) + 3*fop(:, 3);
    
    fop1 = FrenteOptimo(fpareto1, fpareto2, 1, 1);
    
    if fop1(1) > 0
        fop = fop(fop1,:);
    end
    
    %% Fitness y proporcion del mismo
    %fitness 1 max
    fit1 = f1 - (r(1) * g1.^2) - (r(2) * g2.^2);
    neg3 = find(fit1<=0);
    fit1(neg3) = 0;
    tfit1 = sum(fit1);
    selec1 = fit1./tfit1;
    
    %fitness 2 min
    fit2 = f2 + (r(1) * g1.^2) + (r(2) * g2.^2);
    fit2 = 1./(1+fit2);
    tfit2 = sum(fit2);
    selec2 = fit2./tfit2;
    
    %Probabilidad acumulada
    acum1 = cumsum(selec1);
    acum2 = cumsum(selec2);

    %generar numeros aleatorios para conformar la nueva población
    %cruzada por ruleta
    prob1 = rand(100, 1);
    prob2 = rand(100, 1);
    mCromosoma1 = cromosoma;
    mCromosoma2 = cromosoma;
    for i=1:100
        pre1 = find(acum1 >= prob1(i));
        pre2 = find(acum2 >= prob2(i));
        if ~isempty(pre1)
            mCromosoma1(i,:) = cromosoma(pre1(1),:);
        end
        if ~isempty(pre2)
            mCromosoma2(i,:) = cromosoma(pre2(1),:);
        end
    end

    %% Cruce
    % hacer el cruce de cromosomas, el nuevo cromosoma reemplaza al padre #1
    cruce1 = rand(100, 1);
    cruce2 = rand(100, 1);
    padres1 = find(cruce1 < Pc);
    padres2 = find(cruce2 < Pc);
    nCromosoma1 = mCromosoma1;
    nCromosoma2 = mCromosoma2;
    for i=1:length(padres1)-1
        cant = randi(tam-1);
        nCromosoma1(padres1(i),:) = [mCromosoma1(padres1(i), 1:cant) mCromosoma1(padres1(i+1), cant+1:tam)];
    end
    nCromosoma1(length(padres1),:) = [mCromosoma1(padres1(length(padres1)), 1:cant) mCromosoma1(padres1(1), cant+1:tam)];

    for i=1:length(padres2)-1
        cant = randi(tam-1);
        nCromosoma2(padres2(i),:) = [mCromosoma2(padres2(i), 1:cant) mCromosoma2(padres2(i+1), cant+1:tam)];
    end
    nCromosoma2(length(padres2),:) = [mCromosoma2(padres2(length(padres2)), 1:cant) mCromosoma2(padres2(1), cant+1:tam)];
    
    % Calculo de la cantidad de mutaciones
    Tg = 100*tam;
    gM = Tg * Pm;

    %% Mutacion
    % Genes a mutar y su respectiva mutación
    mutar1 = randi(Tg, gM, 1);
    mutar2 = randi(Tg, gM, 1);
    nCromosoma1(mutar1) = ~nCromosoma1(mutar1);
    nCromosoma2(mutar2) = ~nCromosoma2(mutar2);
    
    
    %% evaluacion por separado
    for c=1:100 
       x1(c) = DecimalEscalado(nCromosoma1(c,1:mitad), 0, 32.5);
       x2(c) = DecimalEscalado(nCromosoma2(c,1:mitad), 0, 32.5);
       
       y1(c) = DecimalEscalado(nCromosoma1(c,mitad+1:mitad*2), 0, 43.3);
       y2(c) = DecimalEscalado(nCromosoma2(c,mitad+1:mitad*2), 0, 43.3);
       
       z1(c) = DecimalEscalado(nCromosoma1(c,mitad*2+1:tam), 0, 50);
       z2(c) = DecimalEscalado(nCromosoma2(c,mitad*2+1:tam), 0, 50);
    end
    f_1 = [nCromosoma1 10*x1 + 9*y1 + 8*z1];
    f_1 = sortrows(f_1, 25);
    
    f_2 = [nCromosoma2 10*x2 + 6*y2 + 3*z2];
    f_2 = sortrows(f_2, 25);
    
    cromosoma = [f_2(1:50, 1:tam); f_1(26:75, 1:tam)];
    
    for c=1:100 
       x(c) = DecimalEscalado(cromosoma(c,1:mitad), 0, 32.5);
       y(c) = DecimalEscalado(cromosoma(c,mitad+1:mitad*2), 0, 43.3);
       z(c) = DecimalEscalado(cromosoma(c,mitad*2+1:tam), 0, 50);
    end
    
    %% PROBLEMA A SOLUCIONAR
    %Funcion objetivo 1 a maximizar
    f1 = 10*x + 9*y + 8*z;
    %Funcion objetivo 2 a minimizar
    f2 = 10*x + 6*y + 3*z;
    %Restriccion #1
    g1 = 130 - 4*x - 3*y - 2*z;
    %Restriccion #2
    g2 = 100 - 3*x - 2*y - 2*z;
    
    %% Cierre de iteracion
    desviacion1 = std(f1, 1);
    desviacion2 = std(f2, 1);
end
toc
%% Presentacion de resultados
fpareto1 = 10*fop(:, 1) + 9*fop(:, 2) + 8*fop(:, 3);
fpareto2 = 10*fop(:, 1) + 6*fop(:, 2) + 3*fop(:, 3);
figure
scatter(f1,f2);
hold on
scatter(fpareto1, fpareto2, 'r');
hold off
 