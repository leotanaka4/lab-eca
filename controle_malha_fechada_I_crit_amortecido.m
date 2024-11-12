% Parâmetros do sistema
K = 1.3157;
tau = 0.0894;
K_amp = 5;

syms Ki s
% Função de transferência do motor CC
G = K / (tau*s + 1);
% Função de transferência do controlador PI
C = (Ki/s);

% Função de transferência em malha fechada
L = C * G;
closed_loop_tf = simplify(L / (1 + L));

% Encontrar os coeficientes do polinômio característico do sistema de malha fechada
[~, den] = numden(closed_loop_tf);
den = collect(den, s);  % Coletar os termos em s
coeffs_den = coeffs(den, s, 'All');  % Encontrar os coeficientes

% Normalizar os coeficientes dividindo pelo termo de s^2
coeff_s2 = coeffs_den(1);
coeffs_den = coeffs_den / coeff_s2;

% Fazer delta = 0, ou seja, dois polos no mesmo lugar
a = coeffs_den(1);
b = coeffs_den(2);
c = coeffs_den(3);

delta = b^2 - 4*a*c;

eqns = delta == 0;

% Resolver as equações para Ki
sol = solve(eqns, Ki);

% Valores encontrados para Ki
Ki_val = double(sol); % aq era sol.Ki pq tinha mais de uma variavel
Ki_val = Ki_val/K_amp;

% Função de transferência do controlador PI com os valores encontrados
C_cont = tf([Ki_val], [1, 0]);

% Função de transferência do motor CC
G = tf(K, [tau, 1]);

% Sistema em malha fechada com o controlador PI
closed_loop_cont = feedback(C_cont * K_amp  * G, 1);

step(closed_loop_cont)