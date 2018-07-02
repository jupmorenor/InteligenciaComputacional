%UNIVERSIDAD DISTRITAL FRANCISCO JOSE DE CALDAS
%CIBERNETICA CUALITATIVA 3
%JUAN PABLO MORENO RICO - 20111020059

%OPTIMIZACION POR ENJAMBRE DE PARTICULAS CON RESTRICCIONES
clear

%% Paso 1: poblacion inicial
maxx = 325;
maxy = 433;
maxz = 500;
x(:,1) = randi(maxx,20,1);
y(:,1) = randi(maxy,20,1);
z(:,1) = randi(maxz,20,1);

%% Parametros
a = 150;
alpha = 2;
b = 10;
c = 0.5;
c1 = 1;
c2 = 1;
C = c^alpha;
ncmax = 0.9;
ncmin = 0.4;
imax = 100;
%% Paso 2: evaluar la funcion para la poblacion inicial
f(:,1) = 10*x + 9*y + 8*z;
g1 = 4*x + 3*y + 2*z - 1300;
g2 = 3*x + 2*y + 2*z - 1000;
F = zeros(20,1);
% inicializar la velocidad de las particulas en cero
v = zeros(20, 3);
dv(1) = std(f(:,1));
for i=2:imax
%i=2;
%while (dv > 100)
    %% Paso con restricciones
    q1 = max(0, g1);
    q2 = max(0, g2);
    p1 = a*(1-1./(exp(q1))) + b;
    p2 = a*(1-1./(exp(q2))) + b;
    
    h1 = q1;
    h1(h1>1)=2;
    h1(h1<=1)=1;
    h2 = q2;
    h2(h2>1)=2;
    h2(h2<=1)=1;
    
    C = (c*i)^alpha;
    H = p1.*(q1.^h1) + p2.*(q2.^h2);
    F(:,i-1) = f(:,i-1) - C*H;
    
    %% Paso 4: Hallar el mejor personal por particula
    [max_val, max_pos] = max(F, [], 2);
    pos = sub2ind(size(f), (1:size(f))', max_pos);
    pbest = [x(pos) y(pos) z(pos)];
    % Hallar el mejor global
    max_max = find(F==max(max_val),1);
    gbest = [x(max_max) y(max_max) z(max_max)];
   
    %% Paso 5
    r1 = rand();
    r2 = rand();
    nc = ncmax - (ncmax-ncmin)*i/imax;
    % Calcular las velocidades de cada particula
    par = [x(:,i-1) y(:,i-1) z(:,i-1)];
    v = nc*v + c1*r1*(pbest - par) + c2*r2*(bsxfun(@minus, gbest, par));

    %% Paso 6: Hallar las nuevas posiciones de las particulas
    x(:, i) = x(:, i-1) + v(:,1); 
    y(:, i) = y(:, i-1) + v(:,2); 
    z(:, i) = z(:, i-1) + v(:,3);
    
    % conservando los limites
    x(:, i) = min(max(0, x(:, i)), maxx);
    y(:, i) = min(max(0, y(:, i)), maxy);
    z(:, i) = min(max(0, z(:, i)), maxz);
    % Evaluar la funcion con los nuevos valores
    f(:,i) = 10*x(:,i) + 9*y(:,i) + 8*z(:,i);
    g1 = 4*x(:,i) + 3*y(:,i) + 2*z(:,i) - 1300;
    g2 = 3*x(:,i) + 2*y(:,i) + 2*z(:,i) - 1000;

    %% Paso 7: Revisar convergencia o numero de iteraciones
    dv(i) = std(f(:,i));
    %i = i+1;
end
%scatter3(x(:,i-1), y(:,i-1), z(:,i-1))
figure
hold on
scatter3(x(1,:), y(1,:), z(1,:))
scatter3(x(2,:), y(2,:), z(2,:))
scatter3(x(3,:), y(3,:), z(3,:))
scatter3(x(4,:), y(4,:), z(4,:))
scatter3(x(5,:), y(5,:), z(5,:))
scatter3(x(6,:), y(6,:), z(6,:))
scatter3(x(7,:), y(7,:), z(7,:))
scatter3(x(8,:), y(8,:), z(8,:))
scatter3(x(9,:), y(9,:), z(9,:))
scatter3(x(10,:), y(10,:), z(10,:))
scatter3(x(11,:), y(11,:), z(11,:))
scatter3(x(12,:), y(12,:), z(12,:))
scatter3(x(13,:), y(13,:), z(13,:))
scatter3(x(14,:), y(14,:), z(14,:))
scatter3(x(15,:), y(15,:), z(15,:))
scatter3(x(16,:), y(16,:), z(16,:))
scatter3(x(17,:), y(17,:), z(17,:))
scatter3(x(18,:), y(18,:), z(18,:))
scatter3(x(19,:), y(19,:), z(19,:))
scatter3(x(20,:), y(20,:), z(20,:))
%-----
	
plot3(x(2,:), y(2,:), z(2,:))
plot3(x(3,:), y(3,:), z(3,:))
plot3(x(4,:), y(4,:), z(4,:))
plot3(x(5,:), y(5,:), z(5,:))
plot3(x(6,:), y(6,:), z(6,:))
plot3(x(7,:), y(7,:), z(7,:))
plot3(x(8,:), y(8,:), z(8,:))
plot3(x(9,:), y(9,:), z(9,:))
plot3(x(10,:), y(10,:), z(10,:))
plot3(x(11,:), y(11,:), z(11,:))
plot3(x(12,:), y(12,:), z(12,:))
plot3(x(13,:), y(13,:), z(13,:))
plot3(x(14,:), y(14,:), z(14,:))
plot3(x(15,:), y(15,:), z(15,:))
plot3(x(16,:), y(16,:), z(16,:))
plot3(x(17,:), y(17,:), z(17,:))
plot3(x(18,:), y(18,:), z(18,:))
plot3(x(19,:), y(19,:), z(19,:))
plot3(x(20,:), y(20,:), z(20,:))