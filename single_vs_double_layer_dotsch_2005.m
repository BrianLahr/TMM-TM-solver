clear all; clc; close all;

% =========================================================================
% COMPARAÇÃO SINGLE-LAYER vs DOUBLE-LAYER (Dötsch et al.)
% COM BUSCA ROBUSTA DE MODO FUNDAMENTAL (TM0)
% =========================================================================

% --- 1. Parâmetros Físicos ---
lambda_m = 1.3e-6;          
k_0 = 2*pi / lambda_m;

n_cover = 1.0;   % Ar
n_film  = 2.33;  % Bi:YIG
n_subs  = 1.95;  % GGG

% Faraday Rotation
Theta_F_deg_cm = 1152; 
Theta_F_rad_m  = Theta_F_deg_cm * (pi/180) * 100; 
delta_mag = abs((Theta_F_rad_m * lambda_m * n_film) / pi);

fprintf('=== Parâmetros ===\n');
fprintf('Lambda: %.2e m | Delta: %.5f\n', lambda_m, delta_mag);

% --- 2. Matriz de Rotação (Voigt) ---
rotation_angle = pi/2; 
rot_mat = [1, 0, 0; 0, cos(rotation_angle), sin(rotation_angle); 0, -sin(rotation_angle), cos(rotation_angle)];

% Tensores Isotrópicos
P_subs = rot_mat * diag([n_subs^2, n_subs^2, n_subs^2]) * rot_mat';
P_cov  = rot_mat * diag([n_cover^2, n_cover^2, n_cover^2]) * rot_mat';

% --- 3. Configuração da Varredura ---
thickness_vec = linspace(0.1, 2.0, 1000); % Espessura TOTAL (d_total)
nrps_single = zeros(size(thickness_vec));
nrps_double = zeros(size(thickness_vec));

% Razão d1/d_total para o Double Layer
ratio_d1 = 0.4; 

fprintf('\nIniciando Simulação Comparativa Robusta...\n');

% Loop principal de espessuras
for i = 1:length(thickness_vec)
    d_total = thickness_vec(i) * 1e-6;
    
    % =========================================================
    % CENÁRIO 1: SINGLE LAYER (d1=0, d2=d_total)
    % =========================================================
    t_single = [2e-6; 2e-6 + d_total; 2e-6 + d_total]; 
    
    % Definição dos Tensores Single
    % Forward: +delta
    eps_fwd = [n_film^2, -1j*delta_mag, 0; 1j*delta_mag, n_film^2, 0; 0, 0, n_film^2];
    P_film_s = rot_mat * eps_fwd * rot_mat';
    perms_fwd_s = {P_subs; P_film_s; P_cov};
    
    % Backward: -delta
    eps_bwd = [n_film^2, -1j*(-delta_mag), 0; 1j*(-delta_mag), n_film^2, 0; 0, 0, n_film^2];
    P_film_b_s = rot_mat * eps_bwd * rot_mat';
    perms_bwd_s = {P_subs; P_film_b_s; P_cov};
    
    % Calcula
    nrps_single(i) = calcular_nrps_robusto(k_0, t_single, perms_fwd_s, perms_bwd_s, n_subs, n_film);
    
    % =========================================================
    % CENÁRIO 2: DOUBLE LAYER (+theta / -theta)
    % =========================================================
    d1 = ratio_d1 * d_total;
    d2 = (1 - ratio_d1) * d_total;
    t_double = [2e-6; 2e-6 + d1; 2e-6 + d1 + d2; 2e-6 + d_total];
    
    % Definição dos Tensores Double
    % Forward: Camada 1 (+delta), Camada 2 (-delta)
    eps_p = [n_film^2, -1j*delta_mag, 0; 1j*delta_mag, n_film^2, 0; 0, 0, n_film^2];
    eps_n = [n_film^2, -1j*(-delta_mag), 0; 1j*(-delta_mag), n_film^2, 0; 0, 0, n_film^2];
    P_p = rot_mat * eps_p * rot_mat';
    P_n = rot_mat * eps_n * rot_mat';
    
    perms_fwd_d = {P_subs; P_p; P_n; P_cov}; % L1=Pos, L2=Neg
    perms_bwd_d = {P_subs; P_n; P_p; P_cov}; % L1=Neg, L2=Pos (Inverte tudo)
    
    % Calcula
    nrps_double(i) = calcular_nrps_robusto(k_0, t_double, perms_fwd_d, perms_bwd_d, n_subs, n_film);
    
    % Log progresso
    if mod(i, 10) == 0
        fprintf('d_total = %.2f um | Single: %.2f | Double: %.2f\n', ...
            thickness_vec(i), abs(nrps_single(i)), abs(nrps_double(i)));
    end
end

% --- 4. Plotagem ---
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% Plot Single Layer
plot(thickness_vec, abs(nrps_single), 'b--', 'LineWidth', 2, 'DisplayName', 'Single Layer');
hold on;

% Plot Double Layer
plot(thickness_vec, abs(nrps_double), 'r-', 'LineWidth', 2, 'DisplayName', 'Double Layer (d_1/d_{total}=0.4)');

grid on;
xlabel('Espessura Total do Filme d (\mu m)', 'FontSize', 12);
ylabel('NRPS \Delta\beta (cm^{-1})', 'FontSize', 12);
title('Comparação Single vs Double Layer (Modo Fundamental TM0)', 'FontSize', 14);
subtitle(sprintf('n_{film}=%.2f, \\Theta_F=%d deg/cm, \\lambda=1.3\\mum', n_film, Theta_F_deg_cm));
legend('Location', 'best', 'FontSize', 11);
xlim([0.1 2.0]);

% Ajuste dinâmico do eixo Y
max_y = max([max(abs(nrps_single)), max(abs(nrps_double))]);
if isnan(max_y) || max_y == 0, max_y = 1; end
ylim([0 max_y*1.1]);

% =========================================================================
% FUNÇÃO MESTRE ROBUSTA DE CÁLCULO DE NRPS
% =========================================================================
function val_nrps = calcular_nrps_robusto(k_0, t, perms_fwd, perms_bwd, n_subs, n_film)
    
    num_layers = length(t);
    
    % --- PARTE 1: BUSCA ROBUSTA NO FORWARD ---
    
    % Configuração do Scan
    num_points_scan = 1000; 
    n_scan_min = n_subs + 0.005; 
    n_scan_max = n_film - 0.005; 
    n_eff_scan = linspace(n_scan_min, n_scan_max, num_points_scan);
    beta_scan_vec = n_eff_scan * k_0;
    
    residuals = zeros(size(beta_scan_vec));
    for k = 1:num_points_scan
        val = compute_B1_muller(beta_scan_vec(k), num_layers, t, perms_fwd, k_0);
        residuals(k) = abs(val);
    end
    
    % Encontrar TODOS os mínimos locais
    local_min_idxs = find(islocalmin(residuals));
    found_betas = [];
    
    for idx = local_min_idxs
        if residuals(idx) < 10.0 
            % Define vizinhança
            if idx > 1 && idx < num_points_scan
                b1 = beta_scan_vec(idx-1); b2 = beta_scan_vec(idx); b3 = beta_scan_vec(idx+1);
            else
                b2 = beta_scan_vec(idx); b1 = b2*0.999; b3 = b2*1.001;
            end
            
            % Refina com Muller
            f_obj = @(b) compute_B1_muller(b, num_layers, t, perms_fwd, k_0);
            root = muller(f_obj, b1, b2, b3);
            
            % Valida raiz
            if root ~= 0 && abs(f_obj(root)) < 1e-4
                % Filtra duplicatas
                if isempty(found_betas) || min(abs(found_betas - root)) > 1e-1
                    found_betas = [found_betas, root];
                end
            end
        end
    end
    
    % Ordena para pegar o FUNDAMENTAL (Maior Beta Real)
    if ~isempty(found_betas)
        found_betas = sort(real(found_betas), 'descend');
        beta_fwd = found_betas(1); % <--- GARANTIA DO TM0
    else
        beta_fwd = 0; 
    end
    
    % --- PARTE 2: CALCULAR BACKWARD (SEMEADO PELO FORWARD) ---
    if beta_fwd ~= 0
        % Usa o beta encontrado como chute inicial preciso para o backward
        % Isso impede que ele pule para outro modo, pois a diferença delta é minúscula
        b_seed = beta_fwd;
        f_bwd = @(b) compute_B1_muller(b, num_layers, t, perms_bwd, k_0);
        
        % Intervalo de busca extremamente estreito ao redor do modo forward
        beta_bwd = muller(f_bwd, b_seed*0.99999, b_seed, b_seed*1.00001);
        
        val_nrps = real(beta_fwd - beta_bwd) / 100; % cm^-1
    else
        val_nrps = NaN;
    end
end


% =========================================================================
% FUNÇÕES DO SOLVER (MULLER E COMPUTE_B1)
% =========================================================================
function B1_complex = compute_B1_muller(beta, num_layers, t, layer_perms, k_0)
    A = zeros(num_layers,1);
    B = zeros(num_layers,1);
    A(num_layers) = 0; B(num_layers) = 1;
    for cont = 1:num_layers
        i = (num_layers + 1) - cont; 
        if i ~= 1
            if i == num_layers, w = 0; else, w = t(i) - t(i-1); end
            perm = layer_perms{i};
            eps_x = perm(1,1); eps_y = perm(3,3); delta = imag(perm(1,3));
            perm_m1 = layer_perms{i-1};
            eps_x_m1 = perm_m1(1,1); eps_y_m1 = perm_m1(3,3); delta_m1 = imag(perm_m1(1,3));
            gamma = sqrt(beta^2*(eps_y/eps_x) - k_0^2*((-delta^2)/eps_x + eps_y));
            if real(gamma) < 0, gamma = -gamma; end 
            theta = gamma*w;
            a = beta*delta - gamma*eps_x; b = beta*delta + gamma*eps_x;
            Delta_val = (-delta_m1^2 + eps_x_m1*eps_y_m1)/(-delta^2 + eps_x*eps_y);
            gamma_m1 = sqrt(beta^2*(eps_y_m1/eps_x_m1) - k_0^2*((-delta_m1^2)/eps_x_m1 + eps_y_m1));
            if real(gamma_m1) < 0, gamma_m1 = -gamma_m1; end 
            a_m1 = beta*delta_m1 - gamma_m1*eps_x_m1; b_m1 = beta*delta_m1 + gamma_m1*eps_x_m1;
            exp_neg = exp(-theta); exp_pos = exp(theta);
            A(i-1) = A(i)*exp_neg * ((Delta_val*a - b_m1)/(a_m1 - b_m1)) + B(i)*exp_pos  * ((Delta_val*b - b_m1)/(a_m1 - b_m1));
            B(i-1) = A(i)*exp_neg * ((Delta_val*a - a_m1)/(b_m1 - a_m1)) + B(i)*exp_pos  * ((Delta_val*b - a_m1)/(b_m1 - a_m1));
        end
    end
    B1_complex = B(1); 
end
function f_val = muller(f, x0, x1, x2)
    iter_max = 100; f_tol = 1e-9; x_tol = 1e-9;
    y0 = f(x0); y1 = f(x1); y2 = f(x2);
    iter = 0;
    while(iter <= iter_max)
        iter = iter + 1;
        h1 = x1 - x0; h2 = x2 - x1;
        d1 = (y1 - y0) / h1; d2 = (y2 - y1) / h2;
        a = (d2 - d1) / (h2 + h1); b = a*h2 + d2; c = y2;
        if (a~=0)
            D = sqrt(b*b - 4*a*c);
            if abs(b - D) < abs(b + D), E = b + D; else, E = b - D; end
            dx = -2*c / E;
        elseif (b~=0), dx = -c/b; else, dx = 0; end
        x3 = x2 + dx; x0 = x1; y0 = y1; x1 = x2; y1 = y2; x2 = x3; y2 = f(x2);
        if (abs(dx) < x_tol || abs(y2) < f_tol), break; end
    end
    if (abs(y2) < 1e-4), f_val = x2; else, f_val = 0; end
end