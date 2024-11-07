% Define as matrizes
P1 = [
    -1, 0, 0, 0, 1;
    0, 1, 0, 0, 0;
    3, 1, 0, 0, -1
];

A = [
    -2, 0, 0, 0, 0;
    0, -1, 0, 0, 0;
    1, 0, 1, 0, 1;
    1, 0, 0, 0, 0;
    1, 1, 0, 0, -1
];

P1T = [
    0.5, -0.5, 0.5;
    0, 1, 0;
    0, 0, 0;
    0, 0, 0;
    1.5, -0.5, 0.5
];

% Multiplica as matrizes
result = P1 * A * P1T;

% Exibe o resultado
disp('Resultado de P1 * A * P1^T:')
disp(result)
