clc
clear all
close all

% Variables initiales
x0 =[0.8997; 0.8177; 0.7215; 0.6209; 0.5274; 0.4498 ;0.3910; 0.3195; 0.2418; 0.1673; 0.1037; 0.45]; 
t = 0; 
x = zeros(12, 1);  
u = [13.05; 0.45; 17.4];     
flag = 1;          

delta_x = 1e-4 * ones(size(x));
delta_u = 1e-4 * ones(size(u));  

% Calculer la matrice A et B à partir de x et u
A = calcule_matrice_A(x0, u, delta_x);  
B = calcule_matrice_B(x0, u, delta_u); 



disp(A)

C = zeros(2, 11);  % La matrice C doit maintenant avoir 3 lignes et 12 colonnes
C(1, 1) = 1;  % Mesure de x1
C(2, 11) = 1; % Mesure de x11

n = size(A,1);

% Modification de C pour inclure z_F
C = [C, zeros(size(C, 1), 1)];  % Ajouter une colonne de zéros à C pour correspondre à A avec z_F
C(end+1, end) = 1;  % Ajouter une nouvelle ligne qui mesure seulement z_F (dernier état)


disp(A)
disp(C)


% Recalcul de la matrice W (Observabilité)
W = C;
n = size(A, 1); % Taille de A
for k = 1:n-1
    W = [W; C * (A^k)];
end




ControlMatrix = B; 
for k = 1:n-1
    ControlMatrix = [ControlMatrix, A^k * B];
end

% Vérifier si le système est observable
if rank(W) == n
    disp('Le système est observable.');
else
    disp('Le système n''est pas observable.');
end

% Vérifier si le système est contrôlable
if rank(ControlMatrix) == n
    disp('Le système est contrôlable.');
else
    disp('Le système n''est pas contrôlable.');
end

% Calcul des valeurs propres de la matrice A
valeurs_propres_A = eig(A);

disp('Valeurs propres de A :')
disp(valeurs_propres_A)

if all(real(valeurs_propres_A) < 1)
    disp('Le système est stable en boucle ouverte.');
else
    disp('Le système est instable en boucle ouverte.');
end

% Affichage des matrices
disp('Matrice A :')
disp(A)
disp('Matrice B :')
disp(B)
disp('Matrice C :')
disp(C)
disp('Matrice W (Observabilité) :')
disp(W)
disp('Matrice de contrôlabilité :')
disp(ControlMatrix)


y = C*x0;  % Calcul de la sortie initiale
disp('Sortie initiale y :')
disp(y)

%-------------------- Fonctions de calcul des matrices A et B -------------------

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

function [sys, x0] = modele(t,x,u,flag)


%----------------------------CONSTANTES------------------
% H: rétention molaire
H=[40 10 10 10 10 10 10 10 10 10 80];
% f: numéro du plateau de charge
f=6;
% n: nombre de plateaux
n=2*f-1;
% F: charge
F=10;

%-----------------------------PROGRAMME------------------
if abs(flag) == 1
   
    L = u(1);
    ZF = u(2);
    V = u(3);

    % Calcul des dérivées pour les 11 premiers états
    dx(1) = (1/H(1)) * V * ( k(x(2)) - x(1) );
    for j = 2:f-1
        dx(j) = (1/H(j)) * ( L * ( x(j-1) - x(j) ) + V * ( k(x(j+1)) - k(x(j)) ));
    end
    dx(f) = (1/H(f)) * ( F * (ZF - x(f)) + L * (x(f-1) - x(f)) + V * (k(x(f+1)) - k(x(f))) );
    for j = f+1:n-1
        dx(j) = (1/H(j)) * ( (F+L) * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))) );
    end
    dx(n) = (1/H(n)) * ( (F + L) * (x(n-1) - x(n)) + V * (x(n) - k(x(n))) );

    % Ajout de l'état constant z_f (12ème état)
    dx(12) = 0; % z_f est constant

    sys = [dx];

elseif flag == 0
    % Initialisation des états (ajout de z_f)
    sys = [n + 1 0 n + 3 3 0 0]; 
    x0 = [ 0.8997 0.8177 0.7215 0.6209 0.5274 0.4498 0.3910 0.3195 0.2418 0.1673 0.1037 0.45 ];  % z_f ajouté
elseif flag == 3
    % Sorties : on renvoie uniquement les états souhaités
    sys = [x' 0 0];  % Renvoie uniquement x1, xn, et z_f


else 
    sys = [];
end

end