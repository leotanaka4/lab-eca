data1 = readtable('WFM01.csv'); 
% Extrai as colunas de tempo e canal 1
tempo1 = data1.(1);
canal1 = data1.(2);

% Carrega os dados do arquivo CSV para o Canal 2
data2 = readtable('WFM02.csv');
% Extrai as colunas de tempo e canal 2
tempo2 = data2.(1);
canal2 = data2.(2);
%
% Gera o gráfico de ambos os canais na mesma figura
figure;
plot(tempo1, canal1, 'r-', 'LineWidth', 0.5); % Gráfico do Canal 1 em vermelho
hold on; % Mantém o gráfico atual para adicionar o próximo
plot(tempo2, canal2, 'b-', 'LineWidth', 0.5); % Gráfico do Canal 2 em azul
hold off; % Libera o gráfico

% Personalização do gráfico
xlabel('Tempo (s)');
ylabel('Tensão (V)');
title('Gráfico de Tempo vs. Canais 1 e 2 do Osciloscópio');
legend('Canal 1', 'Canal 2'); % Adiciona uma legenda
grid on;
%}
% Ajustando os valores das regiões estacionárias da onda quadrada positiva
canal1_sup = canal1( tempo1 >= 0.4 & tempo1 <= 0.6);%assume que entra em regime permanente a partir de t = 0.4s
canal2_sup = canal2( tempo1 >= 0.4 & tempo1 <= 0.6);

% Ajustando os valores das regiões estacionárias da onda quadrada positiva
canal1_inf = canal1( tempo1 >= -0.1 & tempo1 <= -0.05);
canal2_inf = canal2( tempo1 >= -0.1 & tempo1 <= -0.05);

%Média dos valores máximos das ondas quadradas
media_canal1_sup = sum(canal1_sup)/(length(canal1_sup));
media_canal2_sup = sum(canal2_sup)/(length(canal2_sup));

%Média dos valores mínimos das ondas quadradas
media_canal1_inf = sum(canal1_inf)/(length(canal1_inf));
media_canal2_inf = sum(canal2_inf)/(length(canal2_inf));

%Valor usado para atingir a região linear editar
step_medio_entrada = (media_canal1_sup+media_canal1_inf)/2;
disp('Valor de excitação para atingir a região linear:');
disp(step_medio_entrada);

%Valor do step de excitação
step_excitacao = media_canal1_sup - media_canal1_inf;
disp('Valor de entrada a partir do ponto de operação na região linear');
disp(step_excitacao);

%Valor da resposta de estado estacionário
step_resposta = media_canal2_sup - media_canal2_inf;
disp('Valor de resposta a partir da aplicação do step');
disp(step_resposta);

%Valor do ganho K
ganho_K = step_resposta/step_excitacao;
disp('Ganho K:');
disp(ganho_K);
disp('V/V')

%plotando os degraus ajustados
t = tempo1(tempo1 > -0.05 & tempo1 < 0.6);
canal1_ajuste = canal1(tempo1 > -0.05 & tempo1 < 0.6);
canal2_ajuste = canal2(tempo1 > -0.05 & tempo1 < 0.6);
%{
% Gera o gráfico de ambos os canais ajustados na mesma figura
figure;
plot(t, canal1_ajuste, 'r-', 'LineWidth', 1.5); % Gráfico do Canal 1 em vermelho
hold on; % Mantém o gráfico atual para adicionar o próximo
plot(t, canal2_ajuste, 'b-', 'LineWidth', 1.5); % Gráfico do Canal 2 em azul
hold off; % Libera o gráfico

% Personalização do gráfico
xlabel('Tempo (s)');
ylabel('Tensão (V)');
title('Gráfico de Tempo vs. Canais 1 e 2 do Osciloscópio');
legend('Canal 1', 'Canal 2'); % Adiciona uma legenda
grid on;
%}
%Centralizando os valores ruidosos do sinal em torno de 0
canal1_zero = canal1_ajuste - media_canal1_inf;
canal2_zero = canal2_ajuste - media_canal2_inf;
%
% Gera o gráfico de ambos os canais centralizados em 0 na mesma figura
figure;
plot(t, canal1_zero, 'r-', 'LineWidth', 0.5); % Gráfico do Canal 1 em vermelho
hold on; % Mantém o gráfico atual para adicionar o próximo
plot(t, canal2_zero, 'b-', 'LineWidth', 0.5); % Gráfico do Canal 2 em azul
hold off; % Libera o gráfico

% Personalização do gráfico
xlabel('Tempo (s)');
ylabel('Tensão (V)');
title('Gráfico de Tempo vs. Canais 1 e 2 do Osciloscópio');
legend('Canal 1', 'Canal 2'); % Adiciona uma legenda
grid on;
%}

%Cálculo da constante de tempo a partir do método da área
canal2_area = (media_canal2_sup - media_canal2_inf) - canal2_zero;
canal2_area = canal2_area(t >= 0 & t <= 0.4);
t_area = t(t >= 0 & t <= 0.4);
A0 = trapz(t_area, canal2_area);
tau_area = A0/(media_canal2_sup - media_canal2_inf);
disp('Constante de tempo:');
disp(tau_area);

G = tf(ganho_K,[tau_area 1]);
figure;
bode(G);


ind=[1:100:length(t)];
tnew=t(ind);
vtnew=canal2_zero(ind);
vanew=canal1_zero(ind);
%plot(tnew, [vanew vtnew]);

%método do logaritmo neperiano
y_inf = media_canal2_sup-media_canal2_inf;
tr_index = find(vtnew >= (y_inf), 1, 'first')-7;
tr = tnew(tr_index);

t_nep = tnew(tnew >= 0 & tnew <= tr);
y_t = vtnew(tnew >= 0 & tnew <= tr);

b = log(y_inf ./ (y_inf - y_t));

a = (t_nep' * b) / (t_nep' * t_nep);

% Determinar tau
tau_nep = 1 / a;

% Exibir o resultado
fprintf('O valor da constante de tempo tau_nep é: %.4f\n', tau_nep);

% simulink :)
degrau.time = t;
degrau.signals.values = canal1_zero;
degrau.signals.dimensions = 1;


load_system('labeca_simulink_2021a');
set_param('labeca_simulink_2021a', 'StopTime', '0.6');
simOut = sim('labeca_simulink_2021a');
disp(simOut);
%close_system('labeca_simulink', 0); 

figure;
plot(t, canal1_zero, 'r-', 'LineWidth', 0.5); % Gráfico do Canal 1 em vermelho
hold on; % Mantém o gráfico atual para adicionar o próximo
plot(t, canal2_zero, 'b-', 'LineWidth', 0.5); % Gráfico do Canal 2 em azul
hold on;
plot(simOut.resposta_area.Time, simOut.resposta_area.Data, 'g-', 'LineWidth', 2.5);
hold on; 
plot(simOut.resposta_nep.Time, simOut.resposta_nep.Data, 'm-', 'LineWidth', 2.5);
legend('Canal 1', 'Canal 2','Resposta Área', 'Resposta Nep');
hold off; % Libera o gráfico

% Personalização do gráfico
xlabel('Tempo (s)');
ylabel('Tensão (V)');
title('Gráfico de Validação dos modelos (Tempo(s) vs Volts(V))');
grid on;

%{
figure;
plot(t_nep, y_t, 'r-', 'LineWidth', 1.5); % Gráfico do Canal 1 em vermelho
hold off; % Libera o gráfico


%fazer intervalo de amostragem t de 10^-3: 600 pontos, pegar a cada 100

ind=[1:100:length(t)];tnew=t(ind);vtnew=canal2_zero(ind);vanew=canal1_zero(ind);
plot(tnew, [vanew vtnew]);
% Personalização do gráfico
xlabel('Tempo (s)');
ylabel('Tensão (V)');
title('Gráfico de Tempo vs. Canais 1 e 2 do Osciloscópio');
legend('Canal 1', 'Canal 2'); % Adiciona uma legenda
grid on;
%}