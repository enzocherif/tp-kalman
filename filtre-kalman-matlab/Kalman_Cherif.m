clear all
clc
close all
% warning('off', 'all');


%% Chargement des données du fichier .mat
% Le fichier contient les matrices At, Bt, Ct, et les données de bruit z
data = load('ExamData_CHERIF.mat');
At = data.At; % Matrice de transition d'état
Bt = data.Bt; % Matrice de commande
Ct = data.Ct; % Matrice de sortie
z = data.z;   % Bruit de mesure (1000 échantillons de dimension 2)

%% Calcul de la moyenne et de l'écart-type du bruit z
% La moyenne est calculée pour chaque dimension (x et y) des échantillons
n = size(z, 2); % Nombre d'échantillons
% La moyenne est calculée séparément pour chaque dimension (x et y)
mean_z_manual = sum(z, 2) / n;

% Affichage des résultats
disp("Moyenne de z calculée manuellement (dimension x, y) :");
disp(mean_z_manual');

%% Vérification avec la fonction MATLAB intégrée
% Utilisation de la fonction 'mean' pour comparaison
mean_z_builtin = mean(z, 2);

disp("Moyenne de z calculée avec la fonction MATLAB (dimension x, y) :");
disp(mean_z_builtin');

% Vérification d'égalité
if isequal(round(mean_z_manual, 6), round(mean_z_builtin, 6))
    disp("Les deux méthodes donnent des résultats identiques !");
else
    disp("Attention, il y a une différence entre les deux calculs.");
end

%% Calcul de l'écart-type manuellement
% Écart-type calculé séparément pour chaque dimension
std_z_manual = sqrt(sum((z - mean_z_manual).^2, 2) / (n - 1));

% Affichage des résultats
disp("Écart-type de z calculé manuellement (dimension x, y) :");
disp(std_z_manual');

%% Vérification avec la fonction MATLAB intégrée
% Utilisation de la fonction 'std' pour comparaison
std_z_builtin = std(z, 0, 2);

disp("Écart-type de z calculé avec la fonction MATLAB (dimension x, y) :");
disp(std_z_builtin');

% Vérification d'égalité
if isequal(round(std_z_manual, 6), round(std_z_builtin, 6))
    disp("Les deux méthodes pour l'écart-type donnent des résultats identiques !");
else
    disp("Attention, il y a une différence entre les deux calculs de l'écart-type.");
end

%% Calcul de la matrice de covariance manuellement
% Initialisation de la matrice de covariance
cov_z_manual = zeros(2, 2);

% Boucle pour accumuler les contributions de chaque échantillon
for i = 1:n
    deviation = z(:, i) - mean_z_manual; % Déviation par rapport à la moyenne
    cov_z_manual = cov_z_manual + deviation * deviation'; % Produit extérieur
end
% Normalisation par n-1
cov_z_manual = cov_z_manual / (n - 1);

% Affichage de la matrice de covariance calculée manuellement
disp("Matrice de covariance calculée manuellement :");
disp(cov_z_manual);

%% Vérification avec la fonction MATLAB intégrée
% Utilisation de la fonction MATLAB 'cov' pour comparer
cov_z_builtin = cov(z'); % Transpose car 'cov' attend les données en colonnes

% Affichage de la matrice de covariance avec MATLAB
disp("Matrice de covariance calculée avec la fonction MATLAB :");
disp(cov_z_builtin);

% Vérification d'égalité
if isequal(round(cov_z_manual, 6), round(cov_z_builtin, 6))
    disp("Les deux méthodes pour la covariance donnent des résultats identiques !");
else
    disp("Attention, il y a une différence entre les deux calculs de la covariance.");
end

cov_z = cov_z_manual;
mean_z = mean_z_manual;

%% Tracé de l'ellipse de confiance
% Niveau de confiance (e.g., 95%)
eta = 0.95;

% Facteur d'échelle pour l'ellipse (distance chi2 pour 2 dimensions)
a = sqrt(-2 * log(1 - eta)); 

% Calcul de l'inverse de la matrice de covariance
inv_cov_z = inv(cov_z);

% Générer des points pour un cercle unitaire
theta = linspace(0, 2*pi, 100); % 100 points entre 0 et 2*pi
cercle_unite = [cos(theta); sin(theta)]; % Points d'un cercle unitaire

% Transformation pour obtenir les points de l'ellipse
% (scaling avec la racine carrée de la covariance et translation par la moyenne)
ellipse_points = a * sqrtm(cov_z) * cercle_unite + mean_z;

%% Tracé
figure;
% Scatter-plot des échantillons
scatter(z(1, :), z(2, :), 10, 'b', 'filled'); % Nuage de points en bleu
hold on;

% Tracé de l'ellipse de confiance
plot(ellipse_points(1, :), ellipse_points(2, :), 'r', 'LineWidth', 2); % Ellipse en rouge

% Configuration du graphique
xlabel('Dimension 1 de z');
ylabel('Dimension 2 de z');
title('Ellipse de Confiance et Nuage de Points du Bruit');
legend('Points du bruit', 'Ellipse de confiance');
grid on;
axis equal; % Échelle égale pour une représentation correcte
hold off;

%% Nombre de points à l'intérieur de l'ellipse
% Calcul des distances de Mahalanobis pour chaque point
mahal_dist = sum(((z - mean_z)' / cov_z) .* (z - mean_z)', 2);

% Seuil basé sur la distribution chi2 pour le niveau de confiance
threshold = chi2inv(eta, 2); % Seuil pour 2 dimensions et niveau de confiance eta

% Vérification si chaque point est dans l'ellipse
points_in_ellipse = mahal_dist <= threshold;

% Calcul du nombre de points dans l'ellipse
num_points_in = sum(points_in_ellipse);
proportion_in = num_points_in / n; % Proportion des points dans l'ellipse

% Affichage des résultats
disp("Nombre total de points :");
disp(n);
disp("Nombre de points dans l'ellipse :");
disp(num_points_in);
disp("Proportion de points dans l'ellipse (en %) :");
disp(proportion_in * 100);

% Comparaison avec le niveau de confiance théorique
disp("Niveau de confiance théorique (en %) :");
disp(eta * 100);

%% Paramètres du filtre de Kalman
% Dimensions du système
nb_etats = size(At, 1); % Nombre d'états
nb_sorties = size(Ct, 1); % Nombre de sorties
nb_iterations = 100; % Nombre d'itérations

% État initial et covariance
etat_reel = zeros(nb_etats, 1); % État réel initial
etat_estime = zeros(nb_etats, 1); % État estimé initial
P_initial = eye(nb_etats); % Covariance initiale de l'estimation

% Définition d'un petit paramètre epsilon
epsilon = 1e-6;

% Signal de commande aléatoire
commande = randn(size(Bt, 2), nb_iterations); % Signal de commande (aléatoire)

% Structure pour stocker les résultats
resultats = struct();

%% Boucle sur les 4 combinaisons de bruit (adapté pour les états 1 et 2 uniquement)
for configuration = 1:4
    % Définir les matrices de covariance Q (bruit de processus) et R (bruit de mesure)
    switch configuration
        case 1
            Q_bruit_processus = epsilon * eye(nb_etats); % Faible bruit de processus
            R_bruit_mesure = epsilon * eye(nb_sorties); % Faible bruit de mesure
        case 2
            Q_bruit_processus = epsilon * eye(nb_etats);
            R_bruit_mesure = 1e-3 * eye(nb_sorties); % Bruit de mesure modéré
        case 3
            Q_bruit_processus = 1e-2 * eye(nb_etats); % Bruit de processus modéré
            R_bruit_mesure = epsilon * eye(nb_sorties);
        case 4
            Q_bruit_processus = 1e-2 * eye(nb_etats); % Bruit de processus modéré
            R_bruit_mesure = 1e-5 * eye(nb_sorties); % Bruit de mesure faible
    end

    % Initialisation des trajectoires et variables (pour tous les états)
    historique_etats_reels = zeros(nb_etats, nb_iterations); % États réels
    historique_etats_estimes = zeros(nb_etats, nb_iterations); % États estimés
    P = P_initial; % Covariance initiale de l'estimation
    historique_MSE = zeros(1, nb_iterations); % Erreur quadratique moyenne
    historique_RMSE = zeros(1, nb_iterations); % Racine carrée de l'erreur quadratique moyenne

    % Simulation pour la configuration actuelle
    for t = 1:nb_iterations
        % Génération de bruit
        bruit_processus = sqrt(diag(Q_bruit_processus)) .* randn(nb_etats, 1);
        bruit_mesure = sqrt(diag(R_bruit_mesure)) .* randn(nb_sorties, 1);

        % Mise à jour de l'état réel
        etat_reel = At * etat_reel + Bt * commande(:, t) + bruit_processus;
        historique_etats_reels(:, t) = etat_reel;

        % Génération de la sortie mesurée
        sortie_mesuree = Ct * etat_reel + bruit_mesure;

        % Phase de prédiction du filtre de Kalman
        etat_prevu = At * etat_estime + Bt * commande(:, t);
        P_prevu = At * P * At' + Q_bruit_processus;

        % Phase de correction du filtre de Kalman
        Gain_Kalman = P_prevu * Ct' / (Ct * P_prevu * Ct' + R_bruit_mesure);
        etat_estime = etat_prevu + Gain_Kalman * (sortie_mesuree - Ct * etat_prevu);
        P = (eye(nb_etats) - Gain_Kalman * Ct) * P_prevu;

        % Stockage des résultats
        historique_etats_estimes(:, t) = etat_estime;

        % Calcul de l'erreur quadratique moyenne (États 1 et 2 seulement)
        erreur = etat_reel(1:2) - etat_estime(1:2); % Différence entre réels et estimés pour x_1, x_2
        historique_MSE(t) = mean(erreur.^2); % Moyenne des erreurs au carré
        historique_RMSE(t) = sqrt(historique_MSE(t)); % Racine carrée de la MSE
    end

    % Sauvegarde des résultats pour x_1 et x_2 uniquement
    resultats(configuration).etats_reels = historique_etats_reels(1:2, :);
    resultats(configuration).etats_estimes = historique_etats_estimes(1:2, :);
    resultats(configuration).MSE_moyenne = mean(historique_MSE);
    resultats(configuration).RMSE_moyenne = sqrt(mean(historique_MSE));
    resultats(configuration).historique_MSE = historique_MSE;
    resultats(configuration).historique_RMSE = historique_RMSE;
end

%% Analyse et visualisation des résultats pour les états 1 et 2 uniquement
% Description des configurations
descriptions = {
    'Sans bruit de mesure / Sans bruit de processus';
    'Avec bruit de mesure / Sans bruit de processus';
    'Sans bruit de mesure / Avec bruit de processus';
    'Avec bruit de mesure / Avec bruit de processus'
};

% Trajectoires des états 1 et 2
for configuration = 1:4
    figure;
    etats_reels = resultats(configuration).etats_reels;
    etats_estimes = resultats(configuration).etats_estimes;
    RMSE_moyenne = resultats(configuration).RMSE_moyenne;

    for dim = 1:2
        subplot(2, 1, dim);
        plot(1:nb_iterations, etats_reels(dim, :), 'b', 'LineWidth', 1.5); hold on;
        plot(1:nb_iterations, etats_estimes(dim, :), 'r--', 'LineWidth', 1.5);
        xlabel('Itérations');
        ylabel(['État ' num2str(dim)]);
        legend('État Réel', 'État Estimé');
        title(['État ' num2str(dim) ' - ' descriptions{configuration}]);
        grid on;
    end
    sgtitle(['États 1 et 2 : ' descriptions{configuration} ' (RMSE = ' num2str(RMSE_moyenne) ')']);
end

%% Tracé des MSE et RMSE pour les états 1 et 2
figure;
for configuration = 1:4
    subplot(2, 2, configuration);
    plot(1:nb_iterations, resultats(configuration).historique_MSE, 'b', 'LineWidth', 1.5);
    xlabel('Itérations');
    ylabel('MSE');
    title(['MSE pour ' descriptions{configuration}]);
    grid on;
end
sgtitle('Erreur Quadratique Moyenne (MSE) pour les États 1 et 2');

figure;
for configuration = 1:4
    subplot(2, 2, configuration);
    plot(1:nb_iterations, resultats(configuration).historique_RMSE, 'r', 'LineWidth', 1.5);
    xlabel('Itérations');
    ylabel('RMSE');
    title(['RMSE pour ' descriptions{configuration}]);
    grid on;
end
sgtitle('Erreur Quadratique Moyenne Racine (RMSE) pour les États 1 et 2');
