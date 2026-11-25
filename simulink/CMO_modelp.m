function [sys, x0] = modele(t, x, u, flag)

%---------------------------DESCRIPTIF-------------------
% Cette "S-function" est utilisée dans le bloc Modèle de Simulink
%
% Entrées: u(1)=L
%          u(2)=V
%          u(3)=ZF
%
% Sorties: sys(1)=x1 
%          sys(2)=xn
%
% Le modèle de la colonne de distillation est représenté par un système d'équations différentielles.
% Ce système est intégré, en continu, sous Simulink, à partir de l'état initial X0.

%----------------------------CONSTANTES------------------
% H: rétention molaire pour chaque plateau
H = [40 10 10 10 10 10 10 10 10 10 80];
% f: numéro du plateau de charge
f = 6;
% n: nombre de plateaux
n = 2 * f - 1;
% F: charge
F = 10;

%-----------------------------PROGRAMME------------------
if abs(flag) == 1
   
    L = u(1);  % Débit de reflux
    ZF = u(2);  % Concentration de la charge
    V = u(3);   % Débit de rebouillage

    % Calcul des dérivées pour les 11 premiers états
    dx(1) = (1/H(1)) * V * (k(x(2)) - x(1));
    for j = 2:f-1
        dx(j) = (1/H(j)) * (L * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
    end
    dx(f) = (1/H(f)) * (F * (ZF - x(f)) + L * (x(f-1) - x(f)) + V * (k(x(f+1)) - k(x(f))));
    for j = f+1:n-1
        dx(j) = (1/H(j)) * ((F+L) * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
    end
    dx(n) = (1/H(n)) * ((F + L) * (x(n-1) - x(n)) + V * (x(n) - k(x(n))));

    sys = [dx];  % Retourner les dérivées des états

elseif flag == 0
    % Initialisation des états
    sys = [n 0 n+2 3 0 0]; 
    x0 = [0.8997 0.8177 0.7215 0.6209 0.5274 0.4498 0.3910 0.3195 0.2418 0.1673 0.1037];

elseif flag == 3
    % Sorties : on renvoie uniquement les états souhaités
    sys = [x(1) x' x(n)];

else 
    sys = [];
end

% Fonction auxiliaire k(x)
function y = k(x)
    alpha = 2;  % Paramètre de volatilité relative
    y = (alpha * x) / (1 + (alpha - 1) * x);  % Fonction de fraction molaire
