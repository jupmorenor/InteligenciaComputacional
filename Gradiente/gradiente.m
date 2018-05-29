%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%METODO DE BUSQUEDA DEL GRADIENTE

clear all

syms z g x X0 x1 x2 t

x = [x1; x2];
z = 2*x1*x2 + 2*x2 - x1^2 - 2*x2^2;
X0 = [0; 0];
e = 0.001;
g = gradient(z);
m = subs(g, x, X0);

while(abs(m(1,1)) > e || abs(m(2,1)) > e)
    L = X0 + subs(g, x, X0)*t;
    f = diff(subs(z, x, L));
    X0 = X0 + subs(g, x, X0)*solve(f, t)
    m = subs(g, x, X0);
end
subs(z, x, X0)
