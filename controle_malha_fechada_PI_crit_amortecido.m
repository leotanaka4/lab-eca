% Parâmetros do sistema
K = 1.3157;
tau = 0.0894;
K_amp = 5;

% Especificações de desempenho
tss = 4*tau;

syms Kp Ti s
% Função de transferência do motor CC
G = K / (tau*s + 1);
% Função de transferência do controlador PI
C = Kp * (1 + 1/(Ti * s));

% Função de transferência em malha fechada
L = C * G;
closed_loop_tf = simplify(L / (1 + L));

% Coeficientes do polinômio característico
[~, den] = numden(closed_loop_tf);
den = collect(den, s); % Coletar os termos em s
coeffs_den = coeffs(den, s, 'All'); % Obter os coeficientes

% Normalizar pelo termo de s^2
coeff_s2 = coeffs_den(1);
coeffs_den = coeffs_den / coeff_s2;

% Extração de coeficientes
a = coeffs_den(1); % Coeficiente de s^2
b = coeffs_den(2); % Coeficiente de s
c = coeffs_den(3); % Coeficiente constante

% Determinar wn e zeta a partir dos coeficientes
wn = sqrt(c); 
zeta = b / (2 * wn);

% Equações de projeto
eq1 = tss == 4 / (zeta * wn); % Tempo de acomodação
delta = b^2 - 4 * a * c; % Delta para verificar polos coincidentes
eq2 = delta == 0; 

% Resolver para Kp e Ti
sol = solve([eq1, eq2], [Kp, Ti]);

% Extrair soluções
Kp_val = double(sol.Kp);
Kp_val = Kp_val / K_amp; % Ajuste pelo ganho do amplificador
Ti_val = double(sol.Ti);

% Função de transferência do controlador PI com os valores encontrados
C_cont = tf([Kp_val*Ti_val, Kp_val], [Ti_val, 0]);

% Função de transferência do motor CC
G = tf(K, [tau, 1]);

% Sistema em malha fechada com o controlador PI
closed_loop_cont = feedback(C_cont * K_amp  * G, 1);

step(closed_loop_cont) % aq tem grau de liberdade suficiente ent da pra ajustar o polinomio direito