clear all; clc; close all;
% =========================================================================
% SOLVER COMPLETO ISOLADOR CROSS STRIP TM (Lohmeyer 2001)
% TMT para beta_0 e beta_1 + Modelo de Overlap para IS/LO
% --- CORRIGIDO PARA EFEITO MO NÃO RECÍPROCO ---
% =========================================================================

fprintf('=== Solver Isolador Cross Strip TM (Lohmeyer 2001) ===\n');

% --- Constantes Físicas ---
eps_0 = 8.8541878188e-12;
mu_0 = 4*pi*1e-7;
c_0 = 1/sqrt(mu_0*eps_0);
lambda_center = 1.3e-6; %
k_0_center = 2*pi/lambda_center;

% --- Parâmetros Geométricos (Tabela 1 do Artigo) ---
t_b = 0.656e-6;      % Espessura Camada Inferior (t_b) [m]
t_t = 0.559e-6;      % Espessura Camada Superior (t_t) [m]
H_core_strip = t_b + t_t; % Altura total do núcleo na faixa (seção bimodal)
t_solver = [t_b; H_core_strip]; % Vetor de limites de camada [t_b; t_b+t_t]
num_layers_v = 4; % Substrato, L1, L2, Cover

% --- Parâmetros de Material (Tabela 1 do Artigo) ---
n_cover = 1.0;     
n_b  = 2.157;      % n_b (Inferior)
n_t  = 2.275;      % n_t (Superior)
n_subs  = 1.950;   % n_s (Substrato)

% --- Magneto-óptica (Tabela 1 do Artigo) ---
% Camada Inferior (L1, espessura t_b)
Theta_F_b_deg_cm = 150;     
Theta_F_b_rad_m  = Theta_F_b_deg_cm * (pi/180) * 100;
% Camada Superior (L2, espessura t_t)
Theta_F_t_deg_cm = -1000;    
Theta_F_t_rad_m  = Theta_F_t_deg_cm * (pi/180) * 100;

% --- CORREÇÃO 1: Remover Rotação Desnecessária e Definir Isotropia ---
rotation_angle = 0; % Não é necessário rotacionar o sistema para Voigt/TM
rot_mat = eye(3); 
P_subs = diag([n_subs^2, n_subs^2, n_subs^2]); % Substrato isotrópico
P_cov  = diag([n_cover^2, n_cover^2, n_cover^2]); % Cobertura isotrópica

% =========================================================================
% PARTE 1: CÁLCULO DAS CONSTANTES DE PROPAGAÇÃO (TMT)
% =========================================================================

fprintf('\n--- 1. Análise Modal (TMT) em Lambda = %.3f um ---\n', lambda_center*1e6);

% --- 1.1 MODO FORWARD (Propagação +z) ---
delta_mag_b_fwd = (n_b * lambda_center * Theta_F_b_rad_m) / pi;
delta_mag_t_fwd = (n_t * lambda_center * Theta_F_t_rad_m) / pi;

% --- CORREÇÃO 2: MO TENSOR ---
% Tensor MO TM/Voigt (M||z, Propagação y): Termos girotópicos em xz (1,3) e zx (3,1)
eps_L1_fwd = [n_b^2, 0, -1j*delta_mag_b_fwd; 
              0, n_b^2, 0; 
              1j*delta_mag_b_fwd, 0, n_b^2];
eps_L2_fwd = [n_t^2, 0, -1j*delta_mag_t_fwd; 
              0, n_t^2, 0; 
              1j*delta_mag_t_fwd, 0, n_t^2];
perms_fwd = {P_subs; eps_L1_fwd; eps_L2_fwd; P_cov}; 

fprintf('   Calculando Modos Forward (0 e 1)...\n');
[neff_fwd_all, beta_fwd_all] = find_mode_robust(k_0_center, num_layers_v, t_solver, perms_fwd, n_subs, n_t);
beta_fwd_0 = beta_fwd_all(1); 
beta_fwd_1 = beta_fwd_all(2); 
neff_fwd_0 = real(beta_fwd_0 / k_0_center);
neff_fwd_1 = real(beta_fwd_1 / k_0_center);

% --- 1.2 MODO BACKWARD (Propagação -z) ---
% Inverter o sinal do termo girotrópico
delta_mag_b_bwd = -delta_mag_b_fwd;
delta_mag_t_bwd = -delta_mag_t_fwd;

% CORREÇÃO 2: MO TENSOR (Novamente, com delta_bwd)
eps_L1_bwd = [n_b^2, 0, -1j*delta_mag_b_bwd; 
              0, n_b^2, 0; 
              1j*delta_mag_b_bwd, 0, n_b^2];
eps_L2_bwd = [n_t^2, 0, -1j*delta_mag_t_bwd; 
              0, n_t^2, 0; 
              1j*delta_mag_t_bwd, 0, n_t^2];
perms_bwd = {P_subs; eps_L1_bwd; eps_L2_bwd; P_cov};

fprintf('   Calculando Modos Backward (0 e 1)...\n');
[neff_bwd_all, beta_bwd_all] = find_mode_refined_bimodal(k_0_center, num_layers_v, t_solver, perms_bwd, [neff_fwd_0, neff_fwd_1]);
beta_bwd_0 = beta_bwd_all(1); 
beta_bwd_1 = beta_bwd_all(2); 
neff_bwd_0 = real(beta_bwd_0 / k_0_center);
neff_bwd_1 = real(beta_bwd_1 / k_0_center);

% --- 1.3 Resultados Chave do TMT ---
fprintf('   > N_eff FWD (0/1): %.7f / %.7f\n', neff_fwd_0, neff_fwd_1);
fprintf('   > N_eff BWD (0/1): %.7f / %.7f\n', neff_bwd_0, neff_bwd_1);

NRPS_diff = (real(beta_fwd_0) - real(beta_bwd_0)) - (real(beta_fwd_1) - real(beta_bwd_1));
L_IS = pi / abs(NRPS_diff); % Comprimento Característico de Isolação (Equação 2)

fprintf('   > Delta NRPS (d0 - d1): %.2e rad/m\n', NRPS_diff);
fprintf('   > Comprimento Característico L_IS: %.3f mm\n', L_IS*1e3);

% =========================================================================
% PARTE 2: ANÁLISE DE ISOLAMENTO (MODELO DE OVERLAP)
% =========================================================================

% --- 2. Parâmetros do Modelo de Overlap (do Artigo) ---
% Coeficientes de sobreposição (overlap) ideais: w0 = w1.
w0 = 0.703; 
w1 = 0.703;

% --- 3. Parâmetros de Comprimento e Plotagem ---
% Análise em torno do comprimento característico L_IS
L_min = L_IS * 0.995; 
L_max = L_IS * 1.005; 
N_pts = 1000;
L_pts = linspace(L_min, L_max, N_pts); 

% --- 4. Diferenças de Propagação (para FWD e BWD) ---
Delta_beta_fwd = beta_fwd_0 - beta_fwd_1;
Delta_beta_bwd = beta_bwd_0 - beta_bwd_1;

% --- 5. Transmissão de Potência (Equação 1 do Paper) ---
P_fwd = w0^2 + w1^2 + 2*w0*w1 .* cos(real(Delta_beta_fwd) * L_pts);
P_bwd = w0^2 + w1^2 + 2*w0*w1 .* cos(real(Delta_beta_bwd) * L_pts);

% --- 6. Cálculo de Isolamento (IS) e Perda (LO) ---
IS = 10 * log10(P_fwd ./ P_bwd);
LO = -10 * log10(P_fwd);

% --- 7. Resultados e Plotagem ---
fprintf('\n=== RESULTADOS (MODELO DE OVERLAP BÍMODAL) ===\n');

[IS_max, L_idx] = max(IS);
L_opt = L_pts(L_idx);
LO_opt = LO(L_idx);

fprintf('Comp. Ótimo (L_opt): %.6f mm\n', L_opt*1e3);
fprintf('Isolamento Máximo: %.2f dB\n', IS_max);
fprintf('Perda no L_opt: %.3f dB\n', LO_opt);

% Plotagem (similar à Fig. 3 do artigo)
figure('Color','w', 'Position', [50, 50, 600, 800]);
% Plot P_fwd e P_bwd (Top)
subplot(3, 1, 1);
plot(L_pts*1e3, P_fwd, 'k-', 'LineWidth', 1.5, 'DisplayName', 'P^f');
hold on;
plot(L_pts*1e3, P_bwd, 'k:', 'LineWidth', 1.5, 'DisplayName', 'P^b');
title('Transmissão de Potência P^f, P^b vs. Comprimento L');
ylabel('P'); legend('show'); grid on;
xline(L_opt*1e3, 'r--', 'LineWidth', 1);
% Plot IS (Middle)
subplot(3, 1, 2);
plot(L_pts*1e3, IS, 'k-', 'LineWidth', 1.5);
title('Isolamento (IS) vs. Comprimento L');
ylabel('IS [dB]'); grid on;
yline(20, 'r--', '20 dB (Alvo)'); 
ylim([-10 40]);
xline(L_opt*1e3, 'r--', 'LineWidth', 1);
% Plot LO (Bottom)
subplot(3, 1, 3);
plot(L_pts*1e3, LO, 'k-', 'LineWidth', 1.5);
title('Perda (LO) vs. Comprimento L');
xlabel('L [mm]'); ylabel('LO [dB]'); grid on;
yline(1, 'r--', '1 dB (Limite)'); 
ylim([0 4]);
xline(L_opt*1e3, 'r--', 'LineWidth', 1);
sgtitle('Desempenho do Isolador (Modelo de Overlap Bimodal)');

% =========================================================================
% FUNÇÕES AUXILIARES TMT/Muller (Necessárias para rodar)
% =========================================================================

function [n_roots, beta_roots] = find_mode_robust(k_0, num_layers, t, layer_perms, n_min, n_max)
    % Ajustado para encontrar os DOIS modos mais altos (0 e 1)
    num_scan_points = 500; 
    n_scan = linspace(n_max-0.0001, n_min+0.0001, num_scan_points); 
    beta_scan = n_scan * k_0;
    residuals = zeros(size(beta_scan));
    for k = 1:length(beta_scan)
        residuals(k) = abs(compute_B1_muller(beta_scan(k), num_layers, t, layer_perms, k_0));
    end
    % Encontra mínimos locais
    is_local_min = [false, (residuals(2:end-1) < residuals(1:end-2)) & (residuals(2:end-1) < residuals(3:end)), false];
    idxs = find(is_local_min);
    
    if isempty(idxs), [~, i_min] = min(residuals); idxs = i_min; end
    
    found_modes_n = [];
    for k = 1:length(idxs)
        beta_guess = beta_scan(idxs(k));
        f_obj = @(b) compute_B1_muller(b, num_layers, t, layer_perms, k_0);
        beta_root = muller(f_obj, beta_guess*0.9995, beta_guess, beta_guess*1.0005);
        n_root = real(beta_root) / k_0;
        
        if abs(f_obj(beta_root)) < 1e-3 && n_root < n_max+0.01 && n_root > n_min-0.01
             % Verifica se o modo é novo
             if isempty(found_modes_n) || min(abs(found_modes_n - n_root)) > 1e-5
                found_modes_n = [found_modes_n, n_root];
             end
        end
    end
    
    if length(found_modes_n) >= 2
        found_modes_n = sort(found_modes_n, 'descend');
        n_roots = found_modes_n(1:2); 
        beta_roots = n_roots * k_0;
    else
        error('A estrutura TMT não suporta dois modos (bimodal). Verifique os parâmetros da Tabela 1.');
    end
end

function [n_roots, beta_roots] = find_mode_refined_bimodal(k_0, num_layers, t, layer_perms, n_guesses)
    % Busca refinada para os modos 0 e 1, usando os resultados FWD como chute inicial
    n_roots = zeros(size(n_guesses));
    for i = 1:length(n_guesses)
        beta_guess = n_guesses(i) * k_0;
        f_obj = @(b) compute_B1_muller(b, num_layers, t, layer_perms, k_0);
        beta_best = muller(f_obj, beta_guess*0.99999, beta_guess, beta_guess*1.00001);
        n_roots(i) = real(beta_best) / k_0;
    end
    beta_roots = n_roots * k_0;
end

function B1_complex = compute_B1_muller(beta, num_layers, t, layer_perms, k_0)
    A = zeros(num_layers,1); B = zeros(num_layers,1);
    A(num_layers) = 0; B(num_layers) = 1; 
    for cont = 1:num_layers
        i = (num_layers + 1) - cont; 
        if i ~= 1
            if i == num_layers, w = 0; elseif i==2, w = t(1); else, w = t(i-1) - t(i-2); end
            
            perm = layer_perms{i}; eps_x = perm(1,1); eps_y = perm(3,3); 
            delta = imag(perm(3,1)); % <<< CORREÇÃO 3: Extrai delta de eps_zx (3,1)
            
            perm_m1 = layer_perms{i-1}; eps_x_m1 = perm_m1(1,1); eps_y_m1 = perm_m1(3,3); 
            delta_m1 = imag(perm_m1(3,1)); % <<< CORREÇÃO 3: Extrai delta_m1 de eps_zx
            
            gamma = sqrt(beta^2*(eps_y/eps_x) - k_0^2*((-delta^2)/eps_x + eps_y));
            if real(gamma) < 0, gamma = -gamma; end 
            theta = gamma*w;
            
            % O TMT para MO Voigt tem a dependência em delta nos termos a e b:
            a = beta*delta - gamma*eps_x; 
            b = beta*delta + gamma*eps_x;
            
            Delta_val = (-delta_m1^2 + eps_x_m1*eps_y_m1)/(-delta^2 + eps_x*eps_y);
            gamma_m1 = sqrt(beta^2*(eps_y_m1/eps_x_m1) - k_0^2*((-delta_m1^2)/eps_x_m1 + eps_y_m1));
            if real(gamma_m1) < 0, gamma_m1 = -gamma_m1; end 
            a_m1 = beta*delta_m1 - gamma_m1*eps_x_m1; 
            b_m1 = beta*delta_m1 + gamma_m1*eps_x_m1;
            
            exp_neg = exp(-theta); exp_pos = exp(theta);
            den = a_m1 - b_m1; if abs(den) < 1e-15, den = 1e-15; end
            
            A(i-1) = A(i)*exp_neg * ((Delta_val*a - b_m1)/den) + B(i)*exp_pos  * ((Delta_val*b - b_m1)/den);
            B(i-1) = A(i)*exp_neg * ((Delta_val*a - a_m1)/(b_m1 - a_m1)) + B(i)*exp_pos  * ((Delta_val*b - a_m1)/(b_m1 - a_m1));
        end
    end
    B1_complex = B(1); 
end

function f_val = muller(f, x0, x1, x2)
    iter_max = 50; f_tol = 1e-10; x_tol = 1e-10;
    y0 = f(x0); y1 = f(x1); y2 = f(x2); iter = 0;
    while(iter <= iter_max)
        iter = iter + 1;
        h1 = x1 - x0; h2 = x2 - x1; d1 = (y1 - y0) / h1; d2 = (y2 - y1) / h2;
        a = (d2 - d1) / (h2 + h1); b = a*h2 + d2; c = y2;
        
        if (abs(a)>1e-15), D = sqrt(b*b - 4*a*c); if abs(b - D) < abs(b + D), E = b + D; else, E = b - D; end; dx = -2*c / E;
        elseif (abs(b)>1e-15), dx = -c/b; else, dx = 0; end
        
        x3 = x2 + dx; x0 = x1; y0 = y1; x1 = x2; y1 = y2; x2 = x3; y2 = f(x2);
        if (abs(dx) < x_tol || abs(y2) < f_tol), break; end
    end
    if (abs(y2) < 1e-4), f_val = x2; else, f_val = 0; end
end
