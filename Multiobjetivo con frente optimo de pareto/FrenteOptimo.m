%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%ALGORITMOS GENETICOS, EJERCICIO DE MULTIPLES OBJETIVOS 
%CRITERIO DE PARETO

%% Frente optimo de pareto
function frenteoptimo = FrenteOptimo(f1, f2, g1, g2)
n = length(f1);
m = length(g1);
v = zeros(n, 1);
for j=1:n
    for k=1:n
        if m==n
            if (f1(j) > f1(k) || f2(j) < f2(k)) && (g1(j) >= 0 && g2(j) >= 0)
                v(j) = v(j) + 1;
            end
        else
            if (f1(j) > f1(k) || f2(j) < f2(k))
                v(j) = v(j) + 1;
            end
        end
    end
end
frenteoptimo = find(max(v)==v);