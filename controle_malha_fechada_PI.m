% Parâmetros do sistema
K = 1.3157;
tau = 0.0894;

% Especificações de desempenho
PO =0.05;
zeta = -log(PO) / sqrt(pi^2 + (log(PO))^2);
tss = 4*tau;
wn = 4/(tss*zeta);
desired_poles = [1, 2*zeta*wn, wn^2];

syms Kp Ti s
% Função de transferência do motor CC
G = K / (tau*s + 1);
% Função de transferência do controlador PI
C = Kp * (1 + 1/(Ti * s));

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

% Polinômio desejado do denominador do sistema de malha fechada
desired_poly = [1, 2*zeta*wn, wn^2];

% Igualar os coeficientes do polinômio característico ao polinômio desejado
eqns = coeffs_den(2:3) == desired_poly(2:3);

% Resolver as equações para Kp e Ti
sol = solve(eqns, [Kp, Ti]);

% Valores encontrados para Kp e Ti
Kp_val = double(sol.Kp);
Ti_val = double(sol.Ti);

% Função de transferência do controlador PI com os valores encontrados
C_cont = tf([Kp_val*Ti_val, Kp_val], [Ti_val, 0]);

% Função de transferência do motor CC
G = tf(K, [tau, 1]);

% Sistema em malha fechada com o controlador PI
closed_loop_cont = feedback(C_cont * G, 1);

step(closed_loop_cont) % aq tem grau de liberdade suficiente ent da pra ajustar o polinomio direito