frequenciesv = [0.200136 ; 0.299556 ; 0.499875 ; 0.7998 ; 1.3961 ; 2.30624 ; 3.85904 ; 6.3798 ; 10.724 ; 17.7425];
c1_folder = "./C1/";
%c2_folder = "./C2/";
import_signals(frequenciesv, 10, c1_folder)%, c2_folder)

%Bode
figure;
w = 2*pi*frequenciesv(1:9);
gjwdb = (20*log10([G_mag_0_20 ; G_mag_0_30 ; G_mag_0_50 ; G_mag_0_80 ; G_mag_1_40 ; G_mag_2_31 ; G_mag_3_86 ; G_mag_6_38;G_mag_10_72]));
g_phase = [G_phase_0_20 ; G_phase_0_30 ; G_phase_0_50 ; G_phase_0_80 ; G_phase_1_40 ; G_phase_2_31 ; G_phase_3_86 ; G_phase_6_38;G_phase_10_72];

% Gráfico de magnitude
subplot(2,1,1); % Cria uma divisão 2x1, selecionando a primeira (magnitude)
semilogx(w, gjwdb, '+');
hold on;
grid on;
p = polyfit(w, gjwdb, 3);
w_fit = linspace(0.2,10.72 , 200);
w_fit = w_fit .* (2 * pi);
gjwdb_fit = polyval(p, w_fit);
semilogx(w_fit, gjwdb_fit);

N = 4;
g_lf = mean(gjwdb_fit(w_fit < w(3)));
g_3db = g_lf - 3;
p_hf = polyfit(log10(w_fit(w_fit > w(8) & w_fit < w(9))), gjwdb_fit(w_fit > w(8) & w_fit < w(9)), 1);
assintota_lf = ones(length(w_fit)) * g_lf;
assintota_3db = ones(length(w_fit)) * g_3db;
assintota_hf = polyval(p_hf, log10(w_fit));

% Encontrando a interseção
diff_curve = gjwdb_fit - assintota_3db;  % Diferença entre o polinômio ajustado e a assíntota
diff_curve = diff_curve(1,:);
idx_cross = find(diff(diff_curve > 0));  % Índice onde a diferença muda de sinal

if ~isempty(idx_cross)
    % Interpolação linear para encontrar a interseção no espaço logarítmico
    log_w_int_bode = interp1(diff_curve(idx_cross:idx_cross+1), log10(w_fit(idx_cross:idx_cross+1)), 0); 
    w_int_bode = 10^log_w_int_bode;  % Voltando para o domínio linear
    y_int_bode = interp1(log10(w_fit), gjwdb_fit, log_w_int_bode);  % Interpolação da magnitude

    % Plotando o ponto de interseção
    semilogx(w_int_bode, y_int_bode, 'ro', 'MarkerSize', 3, 'LineWidth', 2);
    disp(['Frequência de interseção: ', num2str(w_int_bode)]);
else
    disp('Não foi possível encontrar a interseção.');
end

semilogx(w_fit, assintota_lf);
semilogx(w_fit, assintota_3db);
semilogx(w_fit, assintota_hf);
title('Magnitude vs Frequência');
xlabel('Frequência (rad/s)');
ylabel('Magnitude (dB)');
hold off;

%Do gráfico, w = 9.56055
%Cálculo do tau experimental
w_tau = w_int_bode;
tau_bode = 1/w_tau;
K_bode = 10^(g_lf/20);
disp(['tau_bode: ', num2str(tau_bode)])
disp(['K_bode: ', num2str(K_bode)])

% Gráfico de fase
subplot(2,1,2); % Seleciona a segunda parte para o gráfico de fase
semilogx(w, g_phase, '+');
hold on;
grid on;
p_phase = polyfit(w, g_phase, 3);
g_phase_fit = polyval(p_phase, w_fit);
semilogx(w_fit, g_phase_fit);

title('Fase vs Frequência');
xlabel('Frequência (rad/s)');
ylabel('Fase (radianos)');

% Modelo validação
figure;
g_bode = tf(K_bode,[tau_bode 1]);
[mag,phase,wout] = bode(g_bode,w_fit);
mag_graph =20*log10(squeeze(mag));
phase_graph = squeeze(phase)*pi/180;%convertendo pra radianos
subplot(2,1,1);  % Cria um gráfico 2x1 e seleciona o primeiro (magnitude)
semilogx(w, gjwdb, 'r+');  % Plotar dados experimentais
hold on;
semilogx(wout, mag_graph, 'b-', 'LineWidth', 2);  % Plotar magnitude do modelo
grid on;
title('Bode Plot: Magnitude');
xlabel('Frequência (rad/s)');
ylabel('Magnitude (dB)');
legend('Dados experimentais', 'Modelo de validação');
hold off;

subplot(2,1,2);  % Seleciona o segundo gráfico (fase)
semilogx(w, g_phase, '+');  % Plotar dados experimentais de fase
hold on;
semilogx(wout, phase_graph, 'r-', 'LineWidth', 2);  % Plotar fase do modelo
grid on;
title('Bode Plot: Fase');
xlabel('Frequência (rad/s)');
ylabel('Fase (radianos)');
legend('Dados experimentais', 'Modelo de validação');
hold off;



function [a0, c, phi] = fourier_series(f, t, N, T)
    % f: vetor com os valores da função
    % t: vetor com os tempos correspondentes
    % N: número de harmônicos a serem determinados
    
    %T = max(t) - min(t);  % Período
    a0 = calc_a0(f, t, T);  % Calcula o a0 usando trapz
    c = zeros(1, N);       % Vetor para armazenar c_n
    phi = zeros(1, N);     % Vetor para armazenar phi_n
    
    for n = 1:N
        an = calc_an(f, t, T, n);  % Calcula a_n
        bn = calc_bn(f, t, T, n);  % Calcula b_n
        [c(n), phi(n)] = calc_zn(an, bn);  % Calcula c_n e phi_n
    end
end

function a0 = calc_a0(f, t, T)
    % Calcula a_0 usando a regra trapezoidal
    a0 = trapz(t, f) / T;
end

function an = calc_an(f, t, T, n)
    % Calcula a_n para um dado n
    an = (2 / T) * trapz(t, f .* cos(2 * pi * n * t / T));
end

function bn = calc_bn(f, t, T, n)
    % Calcula b_n para um dado n
    bn = (2 / T) * trapz(t, f .* sin(2 * pi * n * t / T));
end

function [cn, phin] = calc_zn(an, bn)
    % Calcula c_n e phi_n a partir de a_n e b_n
    zn = bn + 1j * an;
    cn = abs(zn);
    phin = angle(zn);
end

function [G_mag, G_phase] = calculate_G_single_frequency(u_signal, y_signal, t, N, T)
    % u_signal: vetor de sinal de entrada u_k(t)
    % y_signal: vetor de sinal de saída y_k(t)
    % t: vetor de tempo correspondente aos sinais
    % N: número de harmônicos a serem considerados (normalmente 1)

    % Cálculo dos coeficientes de Fourier para o sinal de entrada
    [a0_u, cu, phiu] = fourier_series(u_signal, t, N, T);
    
    % Cálculo dos coeficientes de Fourier para o sinal de saída
    [a0_y, cy, phiy] = fourier_series(y_signal, t, N, T);

    % Calcular a magnitude e a fase da função de transferência
    G_mag = cy(1) / cu(1);  % Razão entre as magnitudes do primeiro harmônico
    G_phase = phiy(1) - phiu(1);  % Diferença de fases do primeiro harmônico
end

function import_signals(frequencies, num_samples, c1_folder)%, c2_folder)
    % frequencies: vetor de frequências (omega_k) para as amostras
    % num_samples: número de amostras para cada sinal
    % c1_folder: caminho da pasta onde os arquivos csv de c1 estão localizados
    % c2_folder: caminho da pasta onde os arquivos csv de c2 estão localizados

    % Obter lista de arquivos em c1_folder e c2_folder
    c_files = dir(fullfile(c1_folder, '*.csv'));  % Arquivos CSV para signals

    for i = 1:length(frequencies)
        freq = frequencies(i);
        T = 1 / freq;  % Período do sinal
        
        % Obter o nome do arquivo de sinal em ordem
        filename = fullfile(c1_folder, c_files(i).name);
        
        % Exibir os nomes dos arquivos sendo processados
        fprintf('Processando %s para frequência %.2f Hz\n', filename, freq);
        
        
        % Importar os sinais com o tempo na primeira coluna
        data = readmatrix(filename, 'NumHeaderLines', 1);  % Ignora o cabeçalho
        
        % Separar as colunas de sinal e tempo
        t_full = data(:, 1);        % Primeira coluna: tempo
        u_signal_full = data(:, 2); % Segunda coluna: sinal u(t)
        y_signal_full = data(:, 3); % Terceira coluna: sinal y(t)
        %{
        figure;
        hold on;
        plot(t_full, u_signal_full, 'g--', 'LineWidth', 1.5);
        plot(t_full, y_signal_full, 'y--', 'LineWidth', 1.5);
        %}
        
        % Encontrar o índice inicial e final com base no período T e no vetor de tempo
        T0 = t_full(1);  % Tempo inicial
        end_time = T0 + T;
        start_index = find(t_full >= T0, 1, 'first');
        end_index = find(t_full <= end_time, 1, 'last');
        
        % Extrair o intervalo de interesse
        t_segment = t_full(start_index:end_index);
        u_segment = u_signal_full(start_index:end_index);
        y_segment = y_signal_full(start_index:end_index);
        
        % Substituir o ponto decimal por um sublinhado para gerar um nome de variável válido
        freq_str = strrep(sprintf('%.2f', freq), '.', '_');
        
        % Armazenar os sinais em variáveis com nomes específicos
        assignin('base', sprintf('u_%s', freq_str), u_segment);  % Armazena na workspace base
        assignin('base', sprintf('y_%s', freq_str), y_segment);
        assignin('base', sprintf('t_%s', freq_str), t_segment);  % Vetor de tempo cortado
        
        % Descobrir Fourier para o sinal de entrada
        [a0_u, c_u, phi_u] = fourier_series(u_segment, t_segment, 1, T);
        [a0_y, c_y, phi_y] = fourier_series(y_segment, t_segment, 1, T);
        
        % Calcular a magnitude e fase da função de transferência
        [G_mag, G_phase] = calculate_G_single_frequency(u_segment, y_segment, t_segment, 1, T);

        % Armazenar resultados de Fourier e G em variáveis dinâmicas
        assignin('base', sprintf('a0_u_%s', freq_str), a0_u);
        assignin('base', sprintf('c_u_%s', freq_str), c_u);
        assignin('base', sprintf('phi_u_%s', freq_str), phi_u);
        
        assignin('base', sprintf('a0_y_%s', freq_str), a0_y);
        assignin('base', sprintf('c_y_%s', freq_str), c_y);
        assignin('base', sprintf('phi_y_%s', freq_str), phi_y);
        
        assignin('base', sprintf('G_mag_%s', freq_str), G_mag);
        assignin('base', sprintf('G_phase_%s', freq_str), G_phase);
        
        % Calcular as expressões a0 + a1*sin(wt) usando o tempo real para alinhar as senoides
        omega = 2 * pi * str2double(strrep(freq_str, '_', '.'));  % Frequência angular
        %a1_u = c_u(1) * sen(phi_u(1));  % Coeficiente a1 para u
        %a1_y = c_y(1) * sen(phi_y(1));  % Coeficiente a1 para y

        % Correção: usar o tempo real para alinhar as senoides
        u_approx = a0_u + c_u(1) * sin(phi_u(1) + omega*t_segment);
        y_approx = a0_y + c_y(1) * sin(phi_y(1) + omega*t_segment);

        % Criar a figura e plotar as 6 curvas
        
        %
        %Plotar a0, a0+a1 sen wt + phi e o sinal real
        figure;
        hold on;
        plot(t_segment, u_segment, 'r', 'LineWidth', 0.5);  % Sinal real u(t)
        plot(t_segment, y_segment, 'b', 'LineWidth', 0.5);  % Sinal real y(t)
        plot(t_segment, u_approx, 'g-', 'LineWidth', 1.5);  % a0_u + a1_u*sin(wt)
        plot(t_segment, y_approx, 'y-', 'LineWidth', 1.5);  % a0_y + a1_y*sin(wt)
        plot(t_segment, a0_y * ones(size(t_segment)), 'y--', 'LineWidth', 1.5);  % a0_y
        plot(t_segment, a0_u * ones(size(t_segment)), 'g--', 'LineWidth', 1.5);  % a0_u


        %Configurações do gráfico
        legend('u(t)',  'y(t)',  'a_0_u + a_1_u sin(\omega t)',  'a_0_y + a_1_y sin(\omega t)', 'a_0_y','a_0_u', ...
            'Location', 'Best');
        title(sprintf('Comparação de a_0, a_0 + a_1 sin(\\omega t) e Sinais Reais para Frequência %.2f Hz', str2double(strrep(freq_str, '_', '.'))));
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        grid on;
        hold off;        
    end
end