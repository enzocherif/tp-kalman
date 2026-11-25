clc
clear all
close all

%% Initialisation des variables
% Point stationnaire pour x et u
x0 = [0.8997; 0.8177; 0.7215; 0.6209; 0.5274; 0.4498; 0.3910; 0.3195; 0.2418; 0.1673; 0.1037];
u0 = [13.05; 0.45; 17.4]; % L, ZF, V (stationnaire)

% Perturbations pour les calculs numériques
delta_x = 1e-4 * ones(size(x0));
delta_u = 1e-4 * ones(size(u0));

%% Question 2 : Calcul numérique des matrices A et B et C
A_num = calcule_matrice_A(x0, u0, delta_x);
B_num = calcule_matrice_B(x0, u0, delta_u);

disp('Matrice A (numérique) :');
disp(A_num);

disp('Matrice B (numérique) :');
disp(B_num);

% Définition de la matrice C
C = zeros(2, length(x0));
C(1, 1) = 1;  % Mesure de x1
C(2, end) = 1; % Mesure de x11

disp('Matrice C :');
disp(C);

%% Question 3 : Calcul formel des matrices Jacobiennes A et B
% Symboles pour les variables d'état et d'entrée
syms x [11 1] real;
syms u [3 1] real;

% Calcul des équations du modèle (à partir de l'équation dynamique)
dx = modele_symbolique(x, u);

% Calcul des Jacobiennes
A_formel = jacobian(dx, x); % Dérivée partielle par rapport à x
B_formel = jacobian(dx, u); % Dérivée partielle par rapport à u

% Évaluation des Jacobiennes au point stationnaire
A_formel_eval = double(subs(A_formel, [x; u], [x0; u0]));
B_formel_eval = double(subs(B_formel, [x; u], [x0; u0]));

disp('Matrice A (formelle) :');
disp(A_formel_eval);

disp('Matrice B (formelle) :');
disp(B_formel_eval);

%% Vérification des résultats
if norm(A_num - A_formel_eval, 'fro') < 1e-6
    disp('Les matrices A numériques et formelles sont identiques.');
else
    disp('Les matrices A numériques et formelles diffèrent.');
end

if norm(B_num - B_formel_eval, 'fro') < 1e-6
    disp('Les matrices B numériques et formelles sont identiques.');
else
    disp('Les matrices B numériques et formelles diffèrent.');
end

%% Question 4: Vérification de l'observabilité avec le critère de Kalman
% Calcul des matrices A et B
A = A_num;
B = B_num;
% Construction de la matrice d'observabilité
n = size(A_num, 1); % Taille de la matrice d'état
W = C; % Matrice d'observabilité initiale
for k = 1:n-1
    W = [W; C * (A_num^k)]; % Ajout des termes C*A^k
end

% Test d'observabilité
disp('Matrice d''observabilité W :');
disp(W);

if rank(W) == n
    disp('Le système est observable.');
else
    disp('Le système n''est pas observable.');
end

%% Question 5: Vérification de la contrôlabilité avec le critère de Kalman
% Construction de la matrice de contrôlabilité
n = size(A, 1); % Taille de la matrice d'état
ControlMatrix = B; % Matrice de contrôlabilité initiale
for k = 1:n-1
    ControlMatrix = [ControlMatrix, A^k * B]; % Ajout des termes A^k * B
end

% Test de contrôlabilité
disp('Matrice de contrôlabilité :');
disp(ControlMatrix);

if rank(ControlMatrix) == n
    disp('Le système est contrôlable.');
else
    disp('Le système n''est pas contrôlable.');
end

% Vérification de la stabilité en boucle ouverte
% Calcul des valeurs propres de la matrice A
valeurs_propres_A = eig(A);

disp('Valeurs propres de A :');
disp(valeurs_propres_A);

% Test de stabilité
if all(real(valeurs_propres_A) < 0)
    disp('Le système est stable en boucle ouverte.');
else
    disp('Le système est instable en boucle ouverte.');
end

%% Question 6:
% Pas d'échantillonnage
Te = 0.1; % Intervalle de discrétisation (en secondes)

% Calcul des matrices discrètes Ad et Bd
Ad = expm(A * Te); % Exponentielle de matrice pour Ad
Bd = (expm(A * Te) - eye(size(A))) * (A \ B); % Bd dérivé de Ad

% Affichage des résultats
disp('Matrice Ad (système discrétisé) :');
disp(Ad);

disp('Matrice Bd (système discrétisé) :');
disp(Bd);


%% Question 7: 
% Vérification de l'observabilité du système discret
% Construction de la matrice d'observabilité pour le système discret
n = size(Ad, 1); % Taille de la matrice d'état
Wd = C; % Matrice d'observabilité initiale
for k = 1:n-1
    Wd = [Wd; C * (Ad^k)]; % Ajout des termes C*Ad^k
end

% Test d'observabilité
disp('Matrice d''observabilité Wd (système discret) :');
disp(Wd);

if rank(Wd) == n
    disp('Le système discret est observable.');
else
    disp('Le système discret n''est pas observable.');
end

% Vérification de la contrôlabilité du système discret
% Construction de la matrice de contrôlabilité pour le système discret
ControlMatrix_d = Bd; % Matrice de contrôlabilité initiale
for k = 1:n-1
    ControlMatrix_d = [ControlMatrix_d, Ad^k * Bd]; % Ajout des termes Ad^k * Bd
end

% Test de contrôlabilité
disp('Matrice de contrôlabilité ControlMatrix_d (système discret) :');
disp(ControlMatrix_d);

if rank(ControlMatrix_d) == n
    disp('Le système discret est contrôlable.');
else
    disp('Le système discret n''est pas contrôlable.');
end

% Vérification de la stabilité en boucle ouverte
% Calcul des valeurs propres de la matrice Ad
valeurs_propres_Ad = eig(Ad);

disp('Valeurs propres de Ad (système discret) :');
disp(valeurs_propres_Ad);

% Test de stabilité
if all(abs(valeurs_propres_Ad) < 1)
    disp('Le système discret est stable en boucle ouverte.');
else
    disp('Le système discret est instable en boucle ouverte.');
end
%% Fonction pour calculer le modèle dynamique de manière symbolique
function dx = modele_symbolique(x, u)
    % Modèle symbolique pour la colonne de distillation
    H = [40 10 10 10 10 10 10 10 10 10 80];
    f = 6;
    n = 2 * f - 1;
    F = 10;

    % Variables d'entrée
    L = u(1);
    ZF = u(2);
    V = u(3);

    % Calcul des dérivées dynamiques
    dx = sym(zeros(n, 1));
    dx(1) = (1 / H(1)) * V * (k(x(2)) - x(1));
    for j = 2:f-1
        dx(j) = (1 / H(j)) * (L * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
    end
    dx(f) = (1 / H(f)) * (F * (ZF - x(f)) + L * (x(f-1) - x(f)) + V * (k(x(f+1)) - k(x(f))));
    for j = f+1:n-1
        dx(j) = (1 / H(j)) * ((F + L) * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
    end
    dx(n) = (1 / H(n)) * ((F + L) * (x(n-1) - x(n)) + V * (x(n) - k(x(n))));

    function y = k(x)
        alpha = 2;
        y = (alpha * x) / (1 + (alpha - 1) * x);
    end
end

%% Fonctions pour le calcul numérique des matrices A et B
function A = calcule_matrice_A(x0, u0, delta_x)
    n = length(x0);
    A = zeros(n);
    f_x0 = modele(0, x0, u0, 1); 

    for j = 1:n
        x_perturbe = x0;
        x_perturbe(j) = x_perturbe(j) + delta_x(j);
        f_x_perturbe = modele(0, x_perturbe, u0, 1);  
        A(:, j) = (f_x_perturbe - f_x0) / delta_x(j);
    end
end

function B = calcule_matrice_B(x0, u0, delta_u)
    m = length(u0);
    n = length(x0);
    B = zeros(n, m);
    f_u0 = modele(0, x0, u0, 1); 

    for j = 1:m
        u_perturbe = u0;
        u_perturbe(j) = u_perturbe(j) + delta_u(j);
        f_u_perturbe = modele(0, x0, u_perturbe, 1); 
        B(:, j) = (f_u_perturbe - f_u0) / delta_u(j);
    end
end

%% Fonction du modèle dynamique
function [sys, x0] = modele(t, x, u, flag)
    H = [40 10 10 10 10 10 10 10 10 10 80];
    f = 6;
    n = 2 * f - 1;
    F = 10;

    if abs(flag) == 1
        L = u(1);
        ZF = u(2);
        V = u(3);

        dx(1) = (1 / H(1)) * V * (k(x(2)) - x(1));
        for j = 2:f-1
            dx(j) = (1 / H(j)) * (L * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
        end
        dx(f) = (1 / H(f)) * (F * (ZF - x(f)) + L * (x(f-1) - x(f)) + V * (k(x(f+1)) - k(x(f))));
        for j = f+1:n-1
            dx(j) = (1 / H(j)) * ((F + L) * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
        end
        dx(n) = (1 / H(n)) * ((F + L) * (x(n-1) - x(n)) + V * (x(n) - k(x(n))));

        sys = dx';
    elseif flag == 0
        sys = [n 0 n+2 3 0 0];
        x0 = [0.8997 0.8177 0.7215 0.6209 0.5274 0.4498 0.3910 0.3195 0.2418 0.1673 0.1037];
    elseif flag == 3
        sys = [x(1) x' x(n)];
    else
        sys = [];
    end

    function y = k(x)
        alpha = 2;
        y = (alpha * x) / (1 + (alpha - 1) * x);
    end
end
