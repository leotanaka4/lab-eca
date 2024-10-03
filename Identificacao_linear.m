% Dados fornecidos
Va = [0, 0.25389, 0.51509, 0.73795, 1.034, 1.2853, 1.5169, 1.8265, 2.1185, 2.2494, ...
    2.5072, 2.825, 2.9996, 3.2805, 3.5183, 3.6466, 3.8122, 4.5664, 5.2013, 5.5958, ...
    6.1087, 6.686, 7.0504, 7.5012, 8.0346, 8.667, 9.0595, 9.6568, 10.104, 10.504, ...
    11.115, 11.558, 12.166, 12.659, 13.117, 13.54, 14.051, 14.668, 15.081];
Vt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.349, 1.7588, 2.6537, 3.5481, 4.1944, ...
    4.2889, 5.5879, 6.0775, 6.6157, 7.3138, 7.78788, 8.615, 9.4409, 10.133, 10.646, ...
    11.536, 11.993, 12.905, 13.544, 14.291, 14.845, 15.572, 16.428, 17.012];

% Grau do polinômio
n = 4; % Ajuste conforme necessário

% Remover o ponto (0,0) dos dados para o ajuste
Va = Va(1:end);
Vt = Vt(1:end);

% Construir a matriz A para o ajuste polinomial
A = zeros(length(Va), n); % Inicializa a matriz A com zeros
for i = 1:n
    A(:, i) = Va'.^(n-i+1); % Preenche cada coluna com Va^n, Va^(n-1), ..., Va^1
end

% Solução do sistema normal para mínimos quadrados
x = (A' * A) \ (A' * Vt'); % Vetor de coeficientes

% Adicionar o coeficiente correspondente ao ponto (0, 0)
x = [x; 0];

% Exibir os coeficientes ajustados
disp('Coeficientes do polinômio ajustado:')
disp(x')

% Gerar uma curva para visualização do ajuste
Va_plot = linspace(min(Va), max(Va), 100);
Vt_fit = polyval([x(1:end-1)' 0], Va_plot);

% Plotar os pontos de dados e o ajuste
figure;
plot(Va, Vt, 'ro', 'DisplayName', 'Dados originais');
hold on;
plot(Va_plot, Vt_fit, 'b-', 'DisplayName', 'Ajuste Polinomial');
xlabel('Va');
ylabel('Vt');
legend show;
title('Ajuste Polinomial dos Dados Usando Mínimos Quadrados');
grid on;

% Passo 3: Calcular a derivada do polinômio ajustado
p = x;
dp = polyder(p); % Derivada do polinômio ajustado
dVt_dVa = polyval(dp, Va); % Avaliação da derivada ao longo de Va_fit

% Plotando a derivada para identificar a região linear
figure;
plot(Va, dVt_dVa, 'k-', 'LineWidth', 1.5);
xlabel('Va (V)');
ylabel('dVt/dVa');
title('Derivada do Ajuste Polinomial');
grid on;

% Identificação da região linear
linear_region_indices = find(abs(diff(dVt_dVa)) < 0.05); % Critério para identificar estabilidade
linear_region_Va = Va(linear_region_indices);

% Exibição da região linear identificada
disp('Região Linear (Va):');
disp([min(linear_region_Va), max(linear_region_Va)]);
