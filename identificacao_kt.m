% Dados fornecidos
Va = [0, 0.25389, 0.51509, 0.73795, 1.034, 1.2853, 1.5169, 1.8265, 2.1185, 2.2494, ...
    2.5072, 2.825, 2.9996, 3.2805, 3.5183, 3.6466, 3.8122, 4.5664, 5.2013, 5.5958, ...
    6.1087, 6.686, 7.0504, 7.5012, 8.0346, 8.667, 9.0595, 9.6568, 10.104, 10.504, ...
    11.115, 11.558, 12.166, 12.659, 13.117, 13.54, 14.051, 14.668, 15.081];
Vt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.349, 1.7588, 2.6537, 3.5481, 4.1944, ...
    4.2889, 5.5879, 6.0775, 6.6157, 7.3138, 7.78788, 8.615, 9.4409, 10.133, 10.646, ...
    11.536, 11.993, 12.905, 13.544, 14.291, 14.845, 15.572, 16.428, 17.012];
rpm = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84.3, 150, 166, 220, 267.8, 279.7, 343.6, ...
    376.6, 433.2, 457.7, 498.0, 516.7, 582.5, 633.7, 667.0, 712.4, 737.2, 799.2, 824.8, ...
    900.7, 936.6, 980.1, 1032, 1073];

% Filtrando os pontos onde Va é maior que 4.6
Va_nonzero = Va(Va > 4.6);
Vt_nonzero = Vt(Va > 4.6);
rpm_nonzero = rpm(Va > 4.6);

% Conversão de rpm para rad/s
rads_nonzero = rpm_nonzero * (2 * pi / 60); % Convertendo de rpm para rad/s

% Criando a figura com dois subplots
figure;

% Primeiro gráfico: Vt vs rpm
subplot(1, 2, 1); % Subplot 1, na posição da esquerda
scatter(rpm_nonzero, Vt_nonzero, 'b', 'filled'); % Plotando os pontos Vt vs rpm
ylabel('Vt (V)');
xlabel('W (rpm)');
hold on;

% Grau do polinômio
n = 1; % Ajuste conforme necessário

% Construir a matriz A para o ajuste polinomial com rpm
A = zeros(length(rpm_nonzero), n + 1); % Inicializa a matriz A com zeros
for i = 1:n + 1
    A(:, i) = rpm_nonzero'.^(n - i + 1); % Preenche cada coluna com rpm^n, rpm^(n-1), ..., rpm^0
end

% Solução do sistema normal para mínimos quadrados com rpm
x = (A' * A) \ (A' * Vt_nonzero'); % Vetor de coeficientes

% Gerando os dados da linha de regressão com rpm
rpm_fit = linspace(0, max(rpm_nonzero), 100);
Vt_fit = polyval(x, rpm_fit);

% Plotando a linha de regressão linear
plot(rpm_fit, Vt_fit, 'k--', 'LineWidth', 1.5); % Linha de regressão

% Adicionando o valor da inclinação no gráfico com rpm
slope = x(1); % Inclinação da linha (coeficiente angular)
text(0.1 * max(rpm_nonzero), 0.9 * max(Vt_nonzero), sprintf('Kt: %.4f V/rpm', slope), 'FontSize', 12, 'Color', 'k');
grid on;
title('Vt vs W (rpm)');

% Segundo gráfico: Vt vs rads (rad/s)
subplot(1, 2, 2); % Subplot 2, na posição da direita
scatter(rads_nonzero, Vt_nonzero, 'r', 'filled'); % Plotando os pontos Vt vs rads
ylabel('Vt (V)');
xlabel('W (rad/s)');
hold on;

% Construir a matriz A para o ajuste polinomial com rads
A_rads = zeros(length(rads_nonzero), n + 1); % Inicializa a matriz A com zeros
for i = 1:n + 1
    A_rads(:, i) = rads_nonzero'.^(n - i + 1); % Preenche cada coluna com rads^n, rads^(n-1), ..., rads^0
end

% Solução do sistema normal para mínimos quadrados com rads
x_rads = (A_rads' * A_rads) \ (A_rads' * Vt_nonzero'); % Vetor de coeficientes

% Gerando os dados da linha de regressão com rads
rads_fit = linspace(0, max(rads_nonzero), 100);
Vt_fit_rads = polyval(x_rads, rads_fit);

% Plotando a linha de regressão linear
plot(rads_fit, Vt_fit_rads, 'k--', 'LineWidth', 1.5); % Linha de regressão

% Adicionando o valor da inclinação no gráfico com rads
slope_rads = x_rads(1); % Inclinação da linha (coeficiente angular)
text(0.1 * max(rads_nonzero), 0.9 * max(Vt_nonzero), sprintf('Kt: %.4f V/(rad/s)', slope_rads), 'FontSize', 12, 'Color', 'k');
grid on;
title('Vt vs W (rad/s)');

% Exibindo os valores de Kt no terminal
Kt_rpm = 1/slope;
Kt_rads = 1/slope_rads;
fprintf('Kt (V/rpm): %.4f\n', slope);
fprintf('Kt (V/rad/s): %.4f\n', slope_rads);
