clear all; clc; close all;

% =========================================================================
% ANIMAÇÃO: SOLVER GUIA RIB (EIM) VS COMPRIMENTO DE ONDA
% =========================================================================

% --- Configuração da Animação ---
filename_gif = 'animacao_rib_TM00.gif';
lambda_start = 0.5e-6;  % 1.50 um
lambda_stop  = 1.50e-6;  % 1.60 um
num_frames   = 50;       % Quantidade de passos na animação
lambda_vec   = linspace(lambda_start, lambda_stop, num_frames);

fprintf('=== Iniciando Animação (Varredura de Lambda) ===\n');

% --- Constantes Físicas ---
eps_0 = 8.8541878188e-12;
mu_0 = 4*pi*1e-7;
c_0 = 1/sqrt(mu_0*eps_0);

% --- Parâmetros Geométricos (Fixos) ---
size_order = 1e-6;           
H_center = 1.8;  % Altura do núcleo [um]
h_side   = 1.0;  % Altura nas laterais [um]
W_rib    = 2.0;  % Largura do Rib [um]

% --- Materiais (Considerados constantes neste range estreito) ---
n_cover = 1.0;    
n_film  = 2.234;  
n_subs  = 2.214;  
delta_mag = 0; 
rotation_angle = pi/2; 
rot_mat = [1, 0, 0; 0, cos(rotation_angle), sin(rotation_angle); 0, -sin(rotation_angle), cos(rotation_angle)];

% Prepara a Figura
fig = figure('Position', [100, 100, 900, 600], 'Color', 'w');

% =========================================================================
% LOOP DE ANIMAÇÃO
% =========================================================================
for frame_idx = 1:length(lambda_vec)
    
    % Atualiza Lambda e vetores de onda
    lambda_0 = lambda_vec(frame_idx);
    freq = c_0 / lambda_0;
    omega = 2*pi*freq;
    k_0 = 2*pi*freq/c_0;
    
    fprintf('Frame %d/%d: Lambda = %.3f nm... ', frame_idx, num_frames, lambda_0*1e9);
    
    % ---------------------------------------------------------------------
    % 1. Resolver Vertical (Centro)
    % ---------------------------------------------------------------------
    t_center = [0; H_center * size_order; H_center * size_order]; 
    num_layers_v = 3;
    P_subs = rot_mat * diag([n_subs^2, n_subs^2, n_subs^2]) * rot_mat';
    eps_film_tensor = [n_film^2, -1j*delta_mag, 0; 1j*delta_mag, n_film^2, 0; 0, 0, n_film^2];
    P_film = rot_mat * eps_film_tensor * rot_mat';
    P_cov  = rot_mat * diag([n_cover^2, n_cover^2, n_cover^2]) * rot_mat';
    perms_v = {P_subs; P_film; P_cov};
    
    [modos_v_centro, betas_v_centro] = find_all_modes_1D(k_0, num_layers_v, t_center, perms_v, n_subs, n_film);
    
    % ---------------------------------------------------------------------
    % 2. Resolver Vertical (Lateral)
    % ---------------------------------------------------------------------
    t_side = [0; h_side * size_order; h_side * size_order];
    [modos_v_lateral, ~] = find_all_modes_1D(k_0, num_layers_v, t_side, perms_v, n_subs, n_film);
    
    % ---------------------------------------------------------------------
    % 3. EIM Horizontal (Focando apenas no modo Fundamental m=0)
    % ---------------------------------------------------------------------
    if isempty(modos_v_centro)
        fprintf('Cutoff (sem modo vertical).\n');
        continue;
    end
    
    % Pega o modo fundamental vertical (índice 1 no Matlab)
    m = 1; 
    n_eff_core_H = real(modos_v_centro(m));
    beta_vert = betas_v_centro(m);
    
    % Define Cladding Horizontal
    if m <= length(modos_v_lateral)
        n_eff_clad_H = real(modos_v_lateral(m));
    else
        n_eff_clad_H = n_subs;
    end
    
    if n_eff_core_H <= n_eff_clad_H
        fprintf('Cutoff (sem confinamento lateral).\n');
        continue;
    end
    
    % Resolve Horizontal
    t_h = [0; W_rib * size_order; W_rib * size_order];
    num_layers_h = 3;
    P_h_clad = diag([n_eff_clad_H^2, n_eff_clad_H^2, n_eff_clad_H^2]);
    P_h_core = diag([n_eff_core_H^2, n_eff_core_H^2, n_eff_core_H^2]);
    perms_h = {P_h_clad; P_h_core; P_h_clad};
    
    [modos_h, betas_h] = find_all_modes_1D(k_0, num_layers_h, t_h, perms_h, n_eff_clad_H, n_eff_core_H);
    
    if isempty(modos_h)
        fprintf('Cutoff (sem modo horizontal).\n');
        continue;
    end
    
    % Pega o modo fundamental horizontal (índice 1)
    p = 1;
    n_eff_final = modos_h(p);
    beta_final = betas_h(p);
    
    fprintf('OK (n_eff = %.5f)\n', real(n_eff_final));
    
    % ---------------------------------------------------------------------
    % 4. Reconstrução de Campos e Plot
    % ---------------------------------------------------------------------
    
    % Vertical (m=0)
    [Av, Bv, gammav, ~, ~] = processar_modos_encontrados(beta_vert, num_layers_v, t_center, perms_v, eps_0, omega, k_0);
    [x_pts, Hy_v, Ex_v, Ez_v] = calculate_fields_with_margin(num_layers_v, t_center, Av, Bv, gammav, perms_v, omega, eps_0, beta_vert, 2.0e-6);
    
    % Horizontal (p=0)
    [Ah, Bh, gammah, ~, ~] = processar_modos_encontrados(beta_final, num_layers_h, t_h, perms_h, eps_0, omega, k_0);
    [y_pts_raw, Env_h, ~, ~] = calculate_fields_with_margin(num_layers_h, t_h, Ah, Bh, gammah, perms_h, omega, eps_0, beta_final, 3.0e-6);
    y_pts = y_pts_raw - (W_rib * size_order / 2);
    
    % Normalização
    Hy_v_norm = abs(Hy_v) / max(abs(Hy_v));
    Env_h_norm = abs(Env_h) / max(abs(Env_h));
    
    % Malha 2D
    [Y_mesh, X_mesh] = meshgrid(y_pts, x_pts);
    Hy_2D = Hy_v_norm(:) * Env_h_norm(:).';
    
    % --- ATUALIZA O GRÁFICO NA FIGURA EXISTENTE ---
    clf(fig); % Limpa a figura para o próximo frame
    
    % Plot Principal (Hy)
    ax = axes(fig);
    unit_str = 'um';
    X_plot = X_mesh / size_order;
    Y_plot = Y_mesh / size_order;
    
    pcolor(Y_plot, X_plot, abs(Hy_2D)); 
    shading interp; 
    colormap('jet'); 
    c = colorbar; 
    c.Label.String = '|Hy| Normalizado';
    clim([0 1]); % FIXA A ESCALA DE CORES para a animação não piscar
    
    hold on;
    % Desenha Rib
    lw = 2; col = 'w';
    line([-W_rib/2 -W_rib/2], [0 H_center*2], 'Color',col,'LineStyle',':','LineWidth',lw);
    line([ W_rib/2  W_rib/2], [0 H_center*2], 'Color',col,'LineStyle',':','LineWidth',lw);
    yline(H_center, 'w-', 'LineWidth', lw);
    yline(h_side, 'w--', 'LineWidth', 1.5);
    
    title(sprintf('Modo Fundamental TM_{00} @ \\lambda = %.3f nm\nn_{eff} = %.5f', lambda_0*1e9, real(n_eff_final)), 'FontSize', 14);
    xlabel(['Largura y [' unit_str ']'], 'FontSize', 12); 
    ylabel(['Espessura x [' unit_str ']'], 'FontSize', 12);
    colormap('hot');
    axis equal tight; 
    ylim([-0.5, H_center+1.5]); 
    xlim([-W_rib*1.5, W_rib*1.5]);
    
    drawnow;
    
    % --- SALVAR FRAME PARA GIF ---
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if frame_idx == 1
        imwrite(imind, cm, filename_gif, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename_gif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
    
end

fprintf('\nAnimação concluída e salva como "%s".\n', filename_gif);


% =========================================================================
% FUNÇÕES AUXILIARES (Mesmas do kernel anterior)
% =========================================================================

function [found_modes, found_betas] = find_all_modes_1D(k_0, num_layers, t, layer_perms, n_min, n_max)
    num_points = 500; % Reduzido levemente para velocidade da animação
    n_scan = linspace(n_min*0.999, n_max*1.001, num_points);
    beta_scan = n_scan * k_0;
    residuals = zeros(size(beta_scan));
    
    for k = 1:num_points
        val = compute_B1_muller(beta_scan(k), num_layers, t, layer_perms, k_0);
        residuals(k) = abs(val);
    end
    
    local_min_idxs = find(islocalmin(residuals));
    found_modes = []; found_betas = [];
    
    for idx = local_min_idxs
        if residuals(idx) < 10.0
            if idx > 1 && idx < num_points
                b1 = beta_scan(idx-1); b2 = beta_scan(idx); b3 = beta_scan(idx+1);
            else
                b2 = beta_scan(idx); b1 = b2*0.999; b3 = b2*1.001;
            end
            f_obj = @(b) compute_B1_muller(b, num_layers, t, layer_perms, k_0);
            beta_root = muller(f_obj, b1, b2, b3);
            if beta_root ~= 0 && abs(f_obj(beta_root)) < 1e-4
                n_opt = beta_root / k_0;
                if real(n_opt) > n_min + 1e-5
                     if isempty(found_modes) || min(abs(real(found_modes) - real(n_opt))) > 1e-4
                        found_modes = [found_modes, n_opt];
                        found_betas = [found_betas, beta_root];
                     end
                end
            end
        end
    end
    if ~isempty(found_modes)
        [~, sort_idx] = sort(real(found_modes), 'descend');
        found_modes = found_modes(sort_idx); found_betas = found_betas(sort_idx);
    end
end

function [x_points, Hy_points, Ex_points, Ez_points] = calculate_fields_with_margin(num_layers, t, A, B, gamma, layer_perms, omega, eps_0, beta, margin)
    x_min = -margin;
    if t(end) == 0, end_pos = 1e-6; else, end_pos = t(end); end
    x_max = end_pos + margin;
    num_points_total = 400; 
    x_points = linspace(x_min, x_max, num_points_total);
    Hy_points = zeros(size(x_points)); Ex_points = zeros(size(x_points)); Ez_points = zeros(size(x_points));
    
    for idx = 1:length(x_points)
        x = x_points(idx);
        if x < 0
            layer_idx = 1; t_ref = 0;
        elseif x >= end_pos
            layer_idx = num_layers; 
             if num_layers > 1, if t(end) > 0, t_ref = t(end-1); else, t_ref = 0; end, else, t_ref = 0; end
        else
            layer_idx = find(x <= t, 1, 'first');
            if isempty(layer_idx), layer_idx = num_layers; end
            if layer_idx == 1, t_ref = 0; else, t_ref = t(layer_idx); end
        end
        A_l = A(layer_idx); B_l = B(layer_idx); g_l = gamma(layer_idx);
        perm = layer_perms{layer_idx}; eps_x = perm(1,1); eps_y = perm(3,3); delta = imag(perm(1,3));
        Hy_val = A_l * exp(g_l*(x - t_ref)) + B_l * exp(-g_l*(x - t_ref));
        dHy_dx = g_l * A_l * exp(g_l*(x - t_ref)) - g_l * B_l * exp(-g_l*(x - t_ref));
        den = (-delta^2 + eps_x*eps_y);
        if abs(den) < 1e-12, Ex_val=0; Ez_val=0; else
            Ex_val = (1/(omega*eps_0*den)) * (beta*eps_y*Hy_val - delta*dHy_dx);
            Ez_val = (1/(omega*eps_0*den)) * (-1j*eps_x*dHy_dx + 1j*delta*beta*Hy_val);
        end
        Hy_points(idx) = Hy_val; Ex_points(idx) = Ex_val; Ez_points(idx) = Ez_val;
    end
end

% --- Funções Matemáticas ---
function B1_complex = compute_B1_muller(beta, num_layers, t, layer_perms, k_0)
    A = zeros(num_layers,1); B = zeros(num_layers,1);
    A(num_layers) = 0; B(num_layers) = 1;
    for cont = 1:num_layers
        i = (num_layers + 1) - cont; 
        if i ~= 1
            if i == num_layers, w = 0; else, w = t(i) - t(i-1); end
            perm = layer_perms{i}; eps_x = perm(1,1); eps_y = perm(3,3); delta = imag(perm(1,3));
            perm_m1 = layer_perms{i-1}; eps_x_m1 = perm_m1(1,1); eps_y_m1 = perm_m1(3,3); delta_m1 = imag(perm_m1(1,3));
            gamma = sqrt(beta^2*(eps_y/eps_x) - k_0^2*((-delta^2)/eps_x + eps_y));
            if real(gamma) < 0, gamma = -gamma; end 
            theta = gamma*w; a = beta*delta - gamma*eps_x; b = beta*delta + gamma*eps_x;
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
    y0 = f(x0); y1 = f(x1); y2 = f(x2); iter = 0;
    while(iter <= iter_max)
        iter = iter + 1;
        h1 = x1 - x0; h2 = x2 - x1; d1 = (y1 - y0) / h1; d2 = (y2 - y1) / h2;
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
function [A, B, gamma, fator_normalizacao, resultado_integral] = processar_modos_encontrados(beta, num_layers, t, layer_perms, eps_0, omega, k_0)
    A = zeros(num_layers, 1); B = zeros(num_layers, 1); gamma = zeros(num_layers, 1);
    A(num_layers) = 0; B(num_layers) = 1;
    for cont = 1:num_layers
        i = (num_layers + 1) - cont;
        if i ~= 1
            if i==num_layers, w=0; else, w=t(i)-t(i-1); end
            perm = layer_perms{i}; eps_x = perm(1,1); eps_y = perm(3,3); delta = imag(perm(1,3));
            perm_m1 = layer_perms{i-1}; eps_x_m1 = perm_m1(1,1); eps_y_m1 = perm_m1(3,3); delta_m1 = imag(perm_m1(1,3));
            gamma(i) = sqrt(beta^2 * (eps_y/eps_x) - k_0^2 * ((-delta^2)/eps_x + eps_y));
            if real(gamma(i)) < 0, gamma(i) = -gamma(i); end
            theta = gamma(i)*w; a = beta*delta - gamma(i)*eps_x; b = beta*delta + gamma(i)*eps_x;
            Delta_val = (-delta_m1^2 + eps_x_m1*eps_y_m1)/(-delta^2 + eps_x*eps_y);
            gamma_m1 = sqrt(beta^2*(eps_y_m1/eps_x_m1) - k_0^2*((-delta_m1^2)/eps_x_m1 + eps_y_m1));
            if real(gamma_m1) < 0, gamma_m1 = -gamma_m1; end
            a_m1 = beta*delta_m1 - gamma_m1*eps_x_m1; b_m1 = beta*delta_m1 + gamma_m1*eps_x_m1;
            exp_neg = exp(-theta); exp_pos = exp(theta);
            A(i-1) = A(i)*exp_neg*((Delta_val*a-b_m1)/(a_m1-b_m1)) + B(i)*exp_pos*((Delta_val*b-b_m1)/(a_m1-b_m1));
            B(i-1) = A(i)*exp_neg*((Delta_val*a-a_m1)/(b_m1-a_m1)) + B(i)*exp_pos*((Delta_val*b-a_m1)/(b_m1-a_m1));
        else
            perm = layer_perms{i}; eps_x = perm(1,1); eps_y = perm(3,3); delta = imag(perm(1,3));
            gamma(i) = sqrt(beta^2 * (eps_y/eps_x) - k_0^2 * ((-delta^2)/eps_x + eps_y));
            if real(gamma(i)) < 0, gamma(i) = -gamma(i); end
        end
    end
    fator_normalizacao = 1; resultado_integral = 1;
end