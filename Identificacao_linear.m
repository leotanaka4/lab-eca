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

% Construir a matriz A para o ajuste polinomial
A = zeros(length(Va), n); % Inicializa a matriz A com zeros
for i = 1:n
    A(:, i) = Va'.^(n-i+1); % Preenche cada coluna com Va^n, Va^(n-1), ..., Va^1
end

% Solução do sistema normal para mínimos quadrados
x = (A' * A) \ (A' * Vt'); % Vetor de coeficientes

% Exibir os coeficientes ajustados
disp('Coeficientes do polinômio ajustado:')
disp(x')

% Gerar uma curva para visualização do ajuste
Va_plot = linspace(min(Va), max(Va), 100);
Vt_fit = polyval([x' 0], Va_plot);

% Plotar os pontos de dados e o ajuste
figure;
scatter(Va, Vt, 'filled', 'MarkerFaceColor', 'k', 'DisplayName', 'Dados coletados');
hold on;
plot(Va_plot, Vt_fit, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Ajuste Polinomial');
xlabel('Va (V)');
ylabel('Vt (V)');
legend show;
grid on;

% Passo 3: Calcular a derivada do polinômio ajustado
p = x;
dp = polyder([p' 0]); % Derivada do polinômio ajustado
dVt_dVa = polyval(dp, Va); % Avaliação da derivada ao longo de Va

% Plotando a derivada para identificar a região linear
figure;
plot(Va, dVt_dVa, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Derivada do Ajuste Polinomial');
xlabel('Va (V)');
ylabel('dVt/dVa (V/V)');
grid on;
hold on;

% Filtrar os dados para a região linear (4.57 <= Va <= 11.27)
Va_linear = Va(Va >= 4.57 & Va <= 13.117);
Vt_linear = Vt(Va >= 4.57 & Va <= 13.117);

% Realizar a regressão linear na região linear
p_linear = polyfit(Va_linear, Vt_linear, 1);

% Exibir a inclinação da reta (coeficiente angular)
disp('Inclinação da reta na região linear (4.57 V <= Va <= 11.558 V):');
disp(p_linear(1));

% Avaliar os pontos das extremidades no polinômio
Va_extremidades = [4.57, 13.117];
dVt_dVa_extremidades = polyval(dp, Va_extremidades);

% Plotar os pontos que representam as extremidades da região linear
plot(Va_extremidades, dVt_dVa_extremidades, ...
     'mo', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Extremidades da Região Linear');

legend show;

% Plotar a reta constante da inclinação sobre a derivada
plot(Va, p_linear(1)*ones(size(Va)), 'g--', 'LineWidth', 1.5, 'DisplayName', 'K_{reg}');

derivada_regiao_linear = dVt_dVa(Va >= 4.57 & Va <= 13.117);
media = sum(derivada_regiao_linear)/length(derivada_regiao_linear);

plot(Va, media*ones(size(Va)), 'b--', 'LineWidth', 1.5, 'DisplayName', 'K_{mean}');

% Plotar a regressão linear na região linear
figure;
hold on;
scatter(Va_linear, Vt_linear, 'k', 'filled', 'DisplayName', 'Dados Coletados na Região Linear');
plot(Va_linear, polyval(p_linear, Va_linear), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Regressão Linear');
text(5.3, 11.5, sprintf('K_{reg}: %.4f V/V', p_linear(1)), 'FontSize', 12, 'Color', 'k');
xlabel('Va (V)');
ylabel('Vt (V)');
legend show;
grid on;
