clear all
clc
close all
% warning('off', 'all');


data = load('ExamData_CHERIF.mat');

At = data.At;
Bt = data.Bt;
Ct = data.Ct;

z = data.z;

%%%%%%% Calcul de la moyenne %%%%%%
n = size(z,2);
moyenne_z = sum(z,2)/n;
ecart_type_z = std(z, 0, 2);


fprintf("La moyenne de z est de : [%f, %f]\n", moyenne_z(1), moyenne_z(2));


mean_z = mean(z, 2);

fprintf("La moyenne de z avec matlab func est de: [%f, %f]\n", mean_z(1), mean_z(2));
fprintf("L ecart type de z avec matlab func est de: [%f, %f]\n", ecart_type_z(1), ecart_type_z(2));
%%%%%%%% Calcul de la covariance %%%%%%%

cov_z = zeros(2,2);

for i = 1:n
    deviation = z(:,i) -  moyenne_z;
    cov_z = cov_z + deviation*deviation';
end

cov_z = cov_z / (n - 1);


cov_z_verif = cov(z');
fprintf("La matrice de covariance avec la fonction MATLAB est :\n");
disp(cov_z_verif);

%%%%%% Ellipse de confiance %%%%%

eta = 0.95;
a = sqrt(-2*log(1-eta));
inv_cov_z = inv(cov_z);

theta = linspace(0, 2*pi, 100);
cercle_unite = [cos(theta); sin(theta)];

ellipse_points = (mean_z + a*inv_cov_z^(-0.5)) * cercle_unite;

scatter(z(1, :), z(2, :), '.'); 
hold on;
plot(ellipse_points(1, :), ellipse_points(2, :), 'r', 'LineWidth', 1.5);
title('Ellipse de confiance basée sur la définition');
xlabel('Dimension 1 de z');
ylabel('Dimension 2 de z');
saveas(gcf, 'ellipse.png'); % Enregistre en format PNG
hold off;

%%%%% Nombre de point à l interieur de l ellipse %%%%%

compteur = 0;

for i = 1:n
    deviation = z(:, i) - mean_z;    
    mahalanobis_distance = sqrt(deviation' * inv_cov_z * deviation);

    if mahalanobis_distance <= a
        compteur = compteur + 1;
    end
end

nombre_point_interieur = (compteur / n) * 100;

fprintf("Proportion des points dans l'ellipse de confiance : %.2f%%\n", nombre_point_interieur);

fprintf("Niveau de confiance attendu : %.2f%%\n", eta * 100);
fprintf("Niveau de confiance obtenu : %.2f%%\n", nombre_point_interieur);


%%%%% Filtre de Kalman discret %%%%%%

n_iterations = 1000;
n = size(At,1);   % 8
m = size(Ct,1);   % 2

x_estime = zeros(n,1);
P = eye(n); 

t = linspace(0, 10, n_iterations);
u = [sin(t) + 0.1 * (1:n_iterations); 0.5 * cos(t) + 0.05 * (1:n_iterations)];

resultat_simu = struct();
rmse_values = zeros(1, 4);
x_reel = zeros(n,1);
x_reel_listes = zeros(n, n_iterations); 

epsilon = 1e-6;  
for combinaison = 1:4
    switch combinaison
        case 1
            Q = eye(n) * epsilon;  
            R = eye(m) * epsilon;  
            utiliser_bruit_mesure = false;
            utiliser_bruit_processus = false;
        case 2
            Q = eye(n) * epsilon;
            R = eye(m) * 1e-3;    
            utiliser_bruit_mesure = true;
            utiliser_bruit_processus = false;
        case 3
            Q = eye(n) * 1e-2;   
            R = eye(m) * epsilon;


            utiliser_bruit_mesure = false;
            utiliser_bruit_processus = true;
        case 4
            Q = eye(n) * 1e-2;
            R = eye(m) * 1e-5;
            utiliser_bruit_mesure = true;
            utiliser_bruit_processus = true;
    end


    x_estime_listes = zeros(n, n_iterations);
    mse_history = zeros(1, n_iterations);
    rmse_history = zeros(1, n_iterations);
    trace_P_history = zeros(1, n_iterations);
    erreur_etat = zeros(n, n_iterations); 

    for k = 1:n_iterations
    
        if utiliser_bruit_processus
            bruit_processus = sqrt(diag(Q)) .* randn(n, 1);
        else
            bruit_processus = zeros(n, 1);
        end
    
        if utiliser_bruit_mesure
            bruit_mesure = sqrt(diag(R)) .* randn(m, 1);
        else
            bruit_mesure = zeros(m, 1);
        end

        x_reel = At * x_reel + Bt * u(:, k) + bruit_processus;
        x_reel_listes(:, k) = x_reel;
    
        x_pred = At * x_estime + Bt * u(:, k);
        P_pred = At * P * At' + Q;

        K = P_pred * Ct' / (Ct * P_pred * Ct' + R);
    
        y = Ct * x_reel + bruit_mesure;
        x_est = x_pred + K * (y - Ct * x_pred);
        P = (eye(n) - K * Ct) * P_pred;
    
        x_estime = x_est;
    
        x_estime_listes(:, k) = x_est;
        trace_P_history(k) = trace(P);
        erreur_etat(:, k) = x_reel - x_est;
        mse_history(k) = mean(erreur_etat(:, k).^2);
        rmse_history(k) = sqrt(mse_history(k));
    end

    resultat_simu(combinaison).x_estime_listes = x_estime_listes;
    resultat_simu(combinaison).x_reel_listes = x_reel_listes;
    resultat_simu(combinaison).Q = Q;
    resultat_simu(combinaison).R = R;
    resultat_simu(combinaison).MSE = mean(mse_history);
    resultat_simu(combinaison).RMSE = sqrt(resultat_simu(combinaison).MSE);
    resultat_simu(combinaison).mse_history = mse_history;
    resultat_simu(combinaison).rmse_history = rmse_history;
    resultat_simu(combinaison).trace_P_history = trace_P_history;
    resultat_simu(combinaison).erreur_etat = erreur_etat;
end



for combinaison = 1:4
    figure;
    x_reel_listes = resultat_simu(combinaison).x_reel_listes;
    x_estime_listes = resultat_simu(combinaison).x_estime_listes;

    for i = 1:n
        subplot(n, 1, i);
        plot(1:n_iterations, x_reel_listes(i, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'État réel');
        hold on;
        plot(1:n_iterations, x_estime_listes(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'État estimé');
        hold off;
        xlabel('Itérations');
        ylabel(['État ' num2str(i)]);
        title(['Combinaison ' num2str(combinaison) ' - État réel vs État estimé (Dim ' num2str(i) ')']);
        legend;
    end
    sgtitle(['État réel vs État estimé pour la combinaison ' num2str(combinaison)]);
end


figure;

for combinaison = 1:4
    rmse_history = resultat_simu(combinaison).rmse_history;
    subplot(2, 2, combinaison);
    plot(1:n_iterations, rmse_history, 'LineWidth', 1.5);
    xlabel('Itération');
    ylabel('RMSE');
    title(['Combinaison ' num2str(combinaison) ' - RMSE au cours du temps']);
    grid on;
end

sgtitle('RMSE au cours du temps pour toutes les combinaisons');

for combinaison = 1:4
    figure;
    x_reel_listes = resultat_simu(combinaison).x_reel_listes;
    x_estime_listes = resultat_simu(combinaison).x_estime_listes;

    for i = 7:8
        subplot(2, 1, i-6);  
        plot(1:n_iterations, x_reel_listes(i, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'État réel');
        hold on;
        plot(1:n_iterations, x_estime_listes(i, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'État estimé');
        hold off;
        xlabel('Itérations');
        ylabel(['État ' num2str(i)]);
        title(['Combinaison ' num2str(combinaison) ' - État réel vs État estimé (Dim ' num2str(i) ')']);
        legend;
    end
    sgtitle(['État réel vs État estimé pour la combinaison ' num2str(combinaison)]);
end