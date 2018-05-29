%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%BUSQUEDA EN UNA DIMENSION

clear all

syms z dz x x1 x2 x3
z = 12*x - 3*x^4 - 2*x^6
dz = diff(z);
e = 0.0001;
x2 = 0;
x3 = 2;
x1 = (x2+x3)/2;

while (e <= (x3-x2)/2) 
    if (subs(dz,x,x1) >= 0)
        x2 = x1;
    else
        x3 = x1;
    end
    x1 = (x2+x3)/2;
end
x = x1
eval(z)
