%% Sucker-Rod Pump Diagnostic - Everitt's Method
% Based on:
% An Improved Finite-Difference Calculation of Downhole Dynamometer Cards for Sucker-Rod Pumps - Everitt, Jennings (1992)
% Implementation for Uniform Strings
%
% Joe Ercolino - July 2011


%% Warning: Clear all variables
clear all; % close all; home

%% Datos del Pozo
E = 30.5e6; % Módulo de Young (psi)
A = pi/4*(0.75)^2; % Area Transversal de las Cabillas (in^2)
L = 2000; % Longitud de la Sarta de Cabilas (ft)
Gc = 32.2; % Facctor de Coversión (lbm*ft/lbf/s^2)
Rho = 490; % Desidad de las Cabillas (lbm/ft^3)
T = 15; % Velocidad de Bombeo (Ciclos por Minuto)
M = 100; % Número de Nodos - 1 (Nodos: 0 a m)
NT = 360; % Número de Valores de Tiempo Calculados en la Bomba (Un Ciclo Completo con el Primer Punto Repetidos)
C = 0.8; % Coeficiente de Amortiguamiento (s^-1)

%% Carta Dinagráfica de Superficie
load case1s

%% Interpolación de la Carta de Superficie a Inervalos Regulares DT
% La Carta debe estar digitalizada en vectores columnas X y F, subiendo y
% luego bajando. Se admiten puntos repetidos pero no retrocesos.

[Xmin, IXmin] = min(X);
X = circshift(X, 1 - IXmin);
F = circshift(F, 1 - IXmin);
[Xmax, IXmax] = max(X);

Xu = X(1:IXmax);
Xd = [X(IXmax: end); X(1)];
Fu = F(1:IXmax);
Fd = [F(IXmax: end); F(1)];

% Carta Dinamométrica en Ascenso
% plot(Xu, Fu, 'b')
% hold on
% Carta Dinamométrica en Descenso
% plot(Xd, Fd, 'r')
% hold off

S = Xmax - Xmin;

% Cargar Cinemática de la Unidad de Bombeo a Intervalos de Tiempo Regulares (Xi)
% Desplazamiento Positivo en Ascenso

load PUnit

% Ciclo Completo  (Ultimo Punto Igual al Primero)
% Xi = Xi(1:end - 1);

% Interpolación de Xi en NT Puntos
Ni = length(Xi);
Xi = interp1(linspace(0,1,Ni + 1), [Xi; Xi(1)], (0:NT - 1)'/NT);

[Ximin, IXimin] = min(Xi);
Xi = circshift(Xi, 1 - IXimin);
[Ximax, IXimax] = max(Xi);

Xi = Xmin + S*(Xi - Ximin)/(Ximax - Ximin);

Xiu = Xi(1:IXimax);
Xid = Xi(IXimax + 1:end);
Fiu = interp1(Xu, Fu, Xiu, 'linear', 'extrap');
Fid = interp1(Xd, Fd, Xid);
Fi = [Fiu; Fid];

% plot(X, F, 'b.')
% hold on
plot(Xi, Fi, 'b-')
% hold off

%% Valores Calculados
% Incremento de Tiempo (s)
DT = 60/T/(NT - 1);

% Incremento de Longitud (ft)
DX = L/(M - 1);

% Velocidad de Propagación de la Onda (ft/s)
V = sqrt(144*E*Gc/Rho);

% Criterio de Estabilidad
Cest = DX/V/DT;
if Cest >= 1
  disp(['Atención. Coeficiente de Estabilidad: ' num2str(Cest) ' >= 1'])
end


% Coeficiente Alfa
Alpha = (DX/DT^2)*Rho*A/144/Gc;

% Coeficientes de la Ecuación de Ondas Discretizada
A1 = Alpha*(1 + C*DT)/(E*A/DX);
A2 = -(Alpha*(2 + C*DT) - 2*E*A/DX)/(E*A/DX);
A3 = Alpha/(E*A/DX);
A4 = -1;

% Número Total de Valores de Tiempo (Válidos en Superficie)
N = 2*(M -1) + NT; % Número Mínimo de Valores
Nrep = ceil(N/NT); % Número de Repeticiones de la Carta de Superfie
N = Nrep*NT; % Múltiplo de NT - 1

% Matriz Nodal de Desplazamientos: (M + 1)*N U(X, T); (ft)
U = NaN(M + 1, N);

% Inicialización de la Primera Fila de U (Superficie). Conversión a ft.
U(1,:) = repmat(-Xi'/12, 1, Nrep);

% Inicialización de la Segunda Fila de U (Segundo Nodo)
U(2,:) = repmat(Fi'*DX/E/A, 1, Nrep) + U(1,:);

%% Lazo Principal en X (Profundidad)
for i = 3:M + 1
  U(i,i-1:N-i+2) = A1*U(i-1,i:N-i+3) + A2*U(i-1,i-1:N-i+2) + A3*U(i-1,i-2:N-i+1) + A4*U(i-2,i-1:N-i+2);
  %plot(-12*U(i,:), E*A/DX*(U(i,:) - U(i - 1,:)))
  %pause
end

Fp = E*A/DX*(U(M+1,:) - U(M,:));
Xp = -12*(U(M+1,:) - max(U(M+1,:)));

% Xp = NaN(1,N);
% for j = M:N-M+1
%   Xp(j) = -12*((1 + C*DT)*U(M, j + 1) - C*DT*U(M, j) + U(M, j - 1) - U(M - 1, j));
% end
% 
% Fp = NaN(1,N);
% for j = M:N-M+1
%  Fp(j) = (E*A/2/DX)*(3*U(M+1,j) - 4*U(M,j) + U(M-1,j));
% end

hold on
% plot(Up, Fpj,'r.')
plot(Xp, Fp,'r')
hold off

xlabel('Displacement (in)')
ylabel('Load (lbs)')
title('Surface and Downhole Dynamometer Cards')
legend('Surface', 'Downhole')

