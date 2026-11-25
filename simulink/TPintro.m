% Paramètres donnés
R = 8;           % Résistance (Ω)
L = 8e-3;        % Inductance (H)
k = 50;          % Raideur du ressort (N/m)
m = 200e-3;      % Masse (kg)
alpha = 5;       % Frottement visqueux (N·s/m)
beta = 5;        % Facteur de force (N/A)
fv = 0;          % 

% Matrice d'état
A = [-R/L,       0,    -alpha/L;
      0,         0,        1;
     beta/m, -k/m,    -fv/m];

B = [1/L; 0; 0];         % Vecteur d'entrée
C = [0, 1, 0];           % Vecteur de sortie
D = 0;                   % Terme direct

% Affichage des matrices A, B et C
disp('Matrice A :');
disp(A);
disp('Vecteur B :');
disp(B);
disp('Vecteur C :');
disp(C);

% Affichage des valeurs propores de la matrice
disp('valeurs propores Matrice A :');
disp(eig(A));


% Conversion de l'état-espace en fonction de transfert
[num, den] = ss2tf(A, B, C, D);

% Affichage de la fonction de transfert
disp('Numérateur de la fonction de transfert :');
disp(num);
disp('Dénominateur de la fonction de transfert :');
disp(den);

% Création de la fonction de transfert
H = tf(num, den);

% Affichage de la fonction de transfert
disp('Fonction de transfert H(p) :');
H
figure(1)
step(H)

% Visualisation des pôles et des zéros
figure;
pzmap(H);
title('Diagramme des pôles et des zéros');

% Vérification de la stabilité du système
poles = roots(den);
if all(real(poles) < 0)
    disp('Le système est stable.');
else
    disp('Le système est instable.');
end


A1=[0 1;-k/m -(beta*alpha)/(m*R)];
B1=[0;beta/(m+R)];
C1=[1 0];
D1=0;
eig(A1)
hold on 
step(A1,B1,C1,D)