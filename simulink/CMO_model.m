function [sys, x0] = modele(t, x, u, flag)
    %---------------------------DESCRIPTIF-------------------
    % Cette "S-function" est utilisée dans le bloc Modèle de Simulink
    %
    % Entrées: 
    %   u(1) = L  (débit de reflux)
    %   u(2) = V  (débit de rebouillage)
    %   u(3) = ZF (concentration de la charge)
    %
    % Sorties: 
    %   sys(1) = x1 (état en haut de la colonne)
    %   sys(2) = xn (état en bas de la colonne)
    %
    % Le modèle de la colonne de distillation est représenté par un système
    % d'équations différentielles. Ce système est intégré en continu sous
    % Simulink, à partir de l'état initial x0.
    
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
        % Système dynamique (calcul des dérivées)
        
        L = u(1);    % débit de reflux
        ZF = u(2);   % concentration de la charge
        V = u(3);    % débit de rebouillage

        % Calcul des dérivées pour les états du système
        dx(1) = (1/H(1)) * V * (k(x(2)) - x(1));  % Pour le premier plateau (x1)
        
        % Boucle pour les états intermédiaires (2 à f-1)
        for j = 2:f-1
            dx(j) = (1/H(j)) * (L * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
        end
        
        % Plateau de charge (f)
        dx(f) = (1/H(f)) * (F * (ZF - x(f)) + L * (x(f-1) - x(f)) + V * (k(x(f+1)) - k(x(f))));
        
        % Boucle pour les plateaux f+1 à n-1
        for j = f+1:n-1
            dx(j) = (1/H(j)) * ((F + L) * (x(j-1) - x(j)) + V * (k(x(j+1)) - k(x(j))));
        end
        
        % Dernier plateau (n)
        dx(n) = (1/H(n)) * ((F + L) * (x(n-1) - x(n)) + V * (x(n) - k(x(n))));
        
        % Ajout de l'état constant z_f (12ème état)
        dx(12) = 0;  % z_f reste constant
        
        sys = [dx];  % Retourner les dérivées d'états (pour la simulation continue)

    elseif flag == 0
        % Initialisation des états (avec z_f inclus)
        sys = [n + 1 0 n + 3 3 0 0];  % Mise à jour pour inclure z_f
        x0 = [0.8997 0.8177 0.7215 0.6209 0.5274 0.4498 0.3910 0.3195 0.2418 0.1673 0.1037 0.45];  % z_f ajouté à x0

    elseif flag == 3
        % Sorties du système : renvoie x1 et xn
        sys = [x(1); x(12); x(n)];

    else
        sys = [];
    end
end

% Fonction auxiliaire k(x)
function y = k(x)
    alpha = 2;  % Paramètre de volatilité relative
    y = (alpha * x) / (1 + (alpha - 1) * x);  % Fonction de fraction molaire
end

