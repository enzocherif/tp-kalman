function [sys, x0] = modele(t, x, u, flag)

% Vérification des arguments d'entrée
if nargin < 4  % Vérifie si le nombre d'arguments est correct
    error('Not enough input arguments.');
end

%---------------------------DESCRIPTIF-------------------
% Cette "S-function" est utilisée dans le bloc Modèle de Simulink
%
% Entrées: u(1)=L
%              u(2)=V
%              u(3)=ZF
%
% Sorties: sys(1)=x1 
%              sys(2)=xn
%
% Le modèle de la colonne de distillation est représenté par un système d'équations différentielles.
% Ce système est intégré, en continu, sous Simulink, à partir de l'état initial X0

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
   
    L = u(1);    % débit de reflux
    ZF = u(2);   % concentration de la charge
    V = u(3);    % débit de rebouillage

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

    % Ajout de l'état constant z_f (12ème état)
    dx(12) = 0; % z_f est constant

    sys = [dx];

elseif flag == 0
    % Initialisation des états (ajout de z_f)
    sys = [n + 1 0 n + 3 3 0 0];
    x0 = [0.8997 0.8177 0.7215 0.6209 0.5274 0.4498 0.3910 0.3195 0.2418 0.1673 0.1037 0.45];  % z_f ajouté
elseif flag == 3
    % Sorties : on renvoie uniquement les états souhaités
    sys = [x(1); x(2:end); x(end-1); x(end)];

else
    sys = [];
end

end

% Fonction auxiliaire k(x)
function y = k(x)
    alpha = 2;
    y = (alpha * x) / (1 + (alpha - 1) * x);
end
