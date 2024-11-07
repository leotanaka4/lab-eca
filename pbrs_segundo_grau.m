close all

Kt = 0.1527;   %de acordo com o codigo identificacao_Kt (em V/(rad/s))
%   revisar Kt
Ka = 8.9673;  %de acordo com o codigo identificacao_Ka (em (rad/s)/V)
K = 1.3446;    %de acordo com o codigo identificacao_linear (em V/V)
%Kg = 1/K;

Kg = Kt/K;
%Kg = 1/Kt*K;
%Kg = 1/Ka;


data = readtable('WV_2O_5.csv'); 
% Extrai as colunas de tempo e canal 1

% Dados
Va = data.(2);  
Ia = data.(3)*20;   %2.8633;  
Vt = data.(4);
t = data.(1);  

%Va = Va-Va(1);  
%Ia = Ia-Ia(1);  
%Vt = Vt-Vt(1);

V_a.time = t;
V_a.signals.values = Va;
V_a.signals.dimensions = 1;

V_t.time = t;
V_t.signals.values = Vt;
V_t.signals.dimensions = 1;

I_a.time = t;
I_a.signals.values = Ia;
I_a.signals.dimensions = 1;

h = t(2)-t(1);     % intervalo de amostragem

% Cálculo de Ue(t) e Um(t)
Ue = Va - (Kg/Kt) * Vt;  % Equação (4.34a)
Um = Ka * Kt * Ia;       % Equação (4.34b)

%Plotando
figure;
subplot(2,1,1);
plot(t, Ue);
xlabel('Tempo (s)');
ylabel('Ue (V)');
title('Resposta de Ue(t)');

subplot(2,1,2);
plot(t, Um);
xlabel('Tempo (s)');
ylabel('Um (V)');
title('Resposta de Um(t)');

figure;
plot(t, Ia);
hold on;
plot(t, Va);
hold on;
plot(t, Vt);
xlabel('Tempo (s)');
ylabel('Ue (V)');
title('Resposta de Ue(t)');
hold off;


% Construção da matriz Me (equação 4.39a)
n = length(Ia) - 1;  % Número de amostras menos 1 para a construção de Me

Me = [Ia(1:n), Ue(1:n)];  % Montar a matriz Me

% Construção da matriz Mm (equação 4.39b)
Mm = [Vt(1:n), Um(1:n)];  % Montar a matriz Mm

% Agora vamos calcular as soluções pelos mínimos quadrados
% Solução para xe (equação 4.41a)
Xe = inv(Me'*Me)*(Me'*Ia(2:end));  % Mínimos quadrados para xe (4.41a)

% Solução para xm (equação 4.41b)
Xm = inv(Mm'*Mm)*(Mm'*Vt(2:end));  % Mínimos quadrados para xm (4.41b)

% Extração de phi_e, gamma_e, phi_m, gamma_m
phi_e = Xe(1);   % Extração de Φe da solução xe
gamma_e = Xe(2); % Extração de Γe da solução xe

phi_m = Xm(1);   % Extração de Φm da solução xm
gamma_m = Xm(2); % Extração de Γm da solução xm


Ra = (1 - phi_e)/gamma_e;
La = -(Ra*h)/log(phi_e);
f = (1-phi_m)/gamma_m;
J = -(f*h)/log(phi_m);

% Exibição dos resultados
disp('Resultados:');
fprintf('Ra = %.4f\n', Ra);
fprintf('La = %.4f\n', La);
fprintf('f = %.4f\n', f);
fprintf('J = %.4f\n', J);

A = [-Ra/La, -Kg/(Kt*La);
      Ka*Kt/J, -f/J];
B = [1/La;
     0];
C = [0 1];

sys = ss(A,B,C,0);

transferf = tf(sys)
polos = pole(transferf)




