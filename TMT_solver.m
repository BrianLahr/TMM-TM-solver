clear all; clc; close all;

% --- Constantes ---
eps_0 = 8.8541878188e-12; % F/m
mu_0 = 4*pi*1e-7; % H/m
c_0 = 1/sqrt(mu_0*eps_0); % m/s

% --- Entrada de Dados ---
fprintf('=== TMT Solver para Guias Multicamadas (Isotrópico/Anisotrópico/Girotrópico) ===\n');
num_layers = input('Entre com o número de camadas (ex: 3): ');
size_order_layers = input('Entre com a ordem de tamanho das camadas (ex: 1e-6 para microns): ');

% Pré-alocação
layer_thickness = zeros(num_layers, 1);
layer_perms = cell(num_layers, 1);
n_max_structure = 0; % Para definir limite superior da busca
n_min_structure = 100; % Para definir limite inferior

% Rotação Fixa (Geometria Voigt: B em Y)
rotation_angle = pi/2; 
rotation_matrix = [1, 0, 0; 
                   0, cos(rotation_angle), sin(rotation_angle);
                   0, -sin(rotation_angle), cos(rotation_angle)];

for cont = 1:num_layers
    % Loop reverso (N até 1)
    i = (num_layers + 1) - cont; 
    
    fprintf('\n--- Camada %d ---\n', i);
    
    % Tratamento de espessura (igual ao original)
    if i == num_layers
         fprintf('(Topo/Cladding - Espessura considerada semi-infinita)\n');
         th_val = 0; 
         temp_th = input('Espessura (pressione Enter para 0 ou digite valor se for caixa finita): ');
         if ~isempty(temp_th), th_val = temp_th; end
    else
         th_val = input(['Entre com a espessura da camada ' num2str(i) ': ']);
    end
    layer_thickness(i) = th_val * size_order_layers;
    
    fprintf('Tipo do material:\n');
    fprintf('1 - Isotrópico (Aceita complexo, ex: 1.45 - 0.001i)\n');
    fprintf('2 - Anisotrópico (nx, ny, nz)\n');
    fprintf('3 - Girotrópico (Magneto-óptico)\n');
    layer_type = input('Opção: ');
    
    delta = 0;
    if layer_type == 1
        % ATENÇÃO: Agora aceita entrada complexa explicitamente
        n_iso = input(['Índice de refração n da camada ' num2str(i) ' (ex: 3.45 ou 3.45-0.01i): ']);
        
        eps_x = n_iso^2; eps_y = n_iso^2; eps_z = n_iso^2;
        
        % Para o scanner, usamos apenas a parte REAL do índice
        if real(n_iso) > n_max_structure, n_max_structure = real(n_iso); end
        if real(n_iso) < n_min_structure, n_min_structure = real(n_iso); end
        
    elseif layer_type == 2 || layer_type == 3
        fprintf('Nota: Você pode usar números complexos para perdas (ex: 1.5-0.01i)\n');
        n_x = input('n_x: '); n_y = input('n_y: '); n_z = input('n_z: ');
        eps_x = n_x^2; eps_y = n_y^2; eps_z = n_z^2;
        
        % Estimativa baseada apenas na parte REAL
        max_n = max(real([n_x n_y n_z]));
        min_n = min(real([n_x n_y n_z]));
        if max_n > n_max_structure, n_max_structure = max_n; end
        if min_n < n_min_structure, n_min_structure = min_n; end
        
        if layer_type == 3
            delta = input('Constante magneto-óptica (delta): ');
        end
    end
    n_scan_min = n_min_structure * 0.99; 
    n_scan_max = n_max_structure * 1.01;
    
    % Construção do Tensor (suporta complexos automaticamente)
    permittivity_tensor = [eps_x, -1j*delta, 0; 
                           1j*delta, eps_y, 0; 
                           0, 0, eps_z]; 
    layer_perms{i} = rotation_matrix * permittivity_tensor * rotation_matrix'; 
end

% Vetor de interfaces (t)
t = cumsum(layer_thickness);

analysis_type = input(['\nTipo de análise:\n' ...
    '1 - Varredura de Modos (Única Frequência)\n' ...
    '2 - Dispersão (Múltiplas Frequências)\n' ...
    'Sua resposta: ']);

if analysis_type == 1
    % ============================================================
    % ANÁLISE DE FREQUÊNCIA ÚNICA (COM PERDAS/GANHO - MULLER)
    % ============================================================
    
    % Entrada de frequência
    size_order_freq = input('Ordem de grandeza da frequência (ex: 1e9, 1e12): ');   
    freq_val = input('Frequência de análise: ');
    freq = size_order_freq * freq_val;
    
    omega = 2*pi*freq;
    k_0 = 2*pi*freq/c_0; 
    
    fprintf('\n--- Iniciando Varredura de Modos (Scan Real + Muller Complexo) ---\n');
    
    % 1. Varredura inicial no eixo REAL para encontrar candidatos
    num_points = 2000; % Mais pontos para não perder modos estreitos
    n_eff_scan = linspace(n_scan_min, n_scan_max, num_points);
    beta_scan = n_eff_scan * k_0;
    
    residuals = zeros(size(beta_scan));
    
    % Nesta etapa, ainda olhamos para abs() só para achar os "vales"
    for k = 1:num_points
        % Nota: compute_B1_complex retorna o valor complexo. Tiramos abs para o plot inicial.
        err_complex = compute_B1_muller(beta_scan(k), num_layers, t, layer_perms, k_0);
        residuals(k) = abs(err_complex);
    end
    
    % Plot do scan inicial (ajuda a visualizar onde estão as raízes)
    figure;
    semilogy(n_eff_scan, residuals, 'b-', 'LineWidth', 1.5);
    xlabel('Índice Efetivo (Parte Real)');
    ylabel('Erro de Contorno |B_{subs}| (log)');
    title('Varredura Inicial (Busca de Raízes)');
    grid on; xlim([n_scan_min, n_scan_max]);
    
    % Encontrar vales (mínimos locais)
    local_min_idxs = find(islocalmin(residuals));
    
    found_modes = [];     % Armazena n_eff (complexo)
    found_betas = [];     % Armazena beta (complexo)
    
    fprintf('Refinando raízes complexas com Método de Muller...\n');
    
    for idx = local_min_idxs
        % Se o resíduo for "razoavelmente" baixo, tentamos convergir
        if residuals(idx) < 10.0
            
            % Definir 3 pontos iniciais para o Muller baseados no scan
            % (x0, x1, x2 próximos do mínimo encontrado)
            if idx > 1 && idx < num_points
                b1 = beta_scan(idx-1);
                b2 = beta_scan(idx);
                b3 = beta_scan(idx+1);
            else
                % Caso de borda (raro), usa um delta pequeno
                b2 = beta_scan(idx);
                b1 = b2 * 0.999;
                b3 = b2 * 1.001;
            end
            
            % Função anônima para o Muller (B_subs complexo)
            f_obj = @(b) compute_B1_muller(b, num_layers, t, layer_perms, k_0);
            
            % Chama o solver de Muller (código que você enviou)
            beta_root = muller(f_obj, b1, b2, b3);
            
            % Se retornou 0, falhou
            if beta_root ~= 0
                % Verificar se a raiz é válida (erro pequeno)
                erro_final = abs(f_obj(beta_root));
                
                if erro_final < 1e-6
                    n_opt = beta_root / k_0;
                    
                    % Filtra duplicatas (baseado na parte real)
                    if isempty(found_modes) || min(abs(real(found_modes) - real(n_opt))) > 1e-4
                        found_modes = [found_modes, n_opt];
                        found_betas = [found_betas, beta_root];
                        
                        % Plota o ponto encontrado no gráfico
                        hold on;
                        plot(real(n_opt), erro_final, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                    end
                end
            end
        end
    end
    
    % Ordenar por n_eff real decrescente
    if ~isempty(found_modes)
        [~, sort_idx] = sort(real(found_modes), 'descend');
        found_modes = found_modes(sort_idx);
        found_betas = found_betas(sort_idx);
    end
    
    % --- Exibir Resultados ---
    fprintf('\n=== Modos Encontrados (Lossy/Complexo) ===\n');
    for i = 1:length(found_modes)
        n_val = found_modes(i);
        beta_val = found_betas(i);
        
        % Cálculo de Perda (dB/cm)
        % Potência decai como exp(-2*alpha*z), onde alpha = imag(beta)
        % Se imag(beta) for negativo, ajustamos o sinal (convenção física)
        alpha = abs(imag(beta_val)); 
        loss_dB_cm = (alpha * 2 * 4.343) / 100; % converte neper/m para dB/cm
        
        fprintf('Modo %d:\n', i);
        fprintf('  n_eff = %.6f %+.6fi\n', real(n_val), imag(n_val));
        fprintf('  Beta  = %.4e %+.4ei\n', real(beta_val), imag(beta_val));
        fprintf('  Perda = %.4f dB/cm\n', loss_dB_cm);
    end
    
    % =============================================
    % CÁLCULO E PLOT DOS CAMPOS
    % =============================================
    fprintf('\n--- Calculando perfis de campo ---\n');
    resultados_modos = struct();
    
    for modo_idx = 1:length(found_modes)
        beta = found_betas(modo_idx);
        
        % Reutiliza suas funções de processamento (elas funcionam com complexos)
        % Note que A e B serão complexos agora
        [A, B, gamma, fator_normalizacao, resultado_integral] = processar_modos_encontrados(beta, num_layers, t, layer_perms, eps_0, omega, k_0);
        
        [x_points, Hy_points, Ex_points, Ez_points] = calculate_fields(...
            num_layers, t, A, B, gamma, layer_perms, omega, eps_0, beta, 1000);
        
        % Plota usando apenas a parte REAL do n_eff no título
        plot_fields(modo_idx, real(found_modes(modo_idx)), x_points, Hy_points, Ex_points, Ez_points, t, num_layers, size_order_layers);

        % === NOVO: CÁLCULO E PLOT 2D (XY) ===
        fprintf('   -> Gerando Heatmaps 2D (XY)...\n');
        
        % Definir uma largura visual para o eixo Y (ex: 2x a espessura total da estrutura)
        total_thickness = t(end);
        visual_width_y = total_thickness * 2; 
        num_points_y = 100; % Resolução em Y
        
        % Expandir campos
        [X_mesh, Y_mesh, Hy_2D] = expand_fields_to_xy(x_points, Hy_points, visual_width_y, num_points_y);
        [~, ~, Ex_2D] = expand_fields_to_xy(x_points, Ex_points, visual_width_y, num_points_y);
        [~, ~, Ez_2D] = expand_fields_to_xy(x_points, Ez_points, visual_width_y, num_points_y);
        
        % Plotar
        plot_2D_heatmaps(modo_idx, real(found_modes(modo_idx)), X_mesh, Y_mesh, Hy_2D, Ex_2D, Ez_2D, t, size_order_layers);
        % ====================================
    end

elseif analysis_type == 2
    % ============================================================
    % ANÁLISE DE DISPERSÃO (MÚLTIPLAS FREQUÊNCIAS)
    % ============================================================
    
    % --- Entrada de dados para análise dispersiva ---
    fprintf('\n=== Configuração da Análise de Dispersão ===\n');
    
    % Solicitar range de frequências
    size_order_freq = input('Ordem de grandeza da frequência (ex: 1e9, 1e12): ');
    fprintf('Defina o range de frequências:\n');
    f_min_val = input('Frequência mínima: ');
    f_max_val = input('Frequência máxima: ');
    num_freq_points = input('Número de pontos de frequência: ');
    
    f_min = f_min_val * size_order_freq;
    f_max = f_max_val * size_order_freq;

    % Criar vetor de frequências
    freq_vec = linspace(f_min, f_max, num_freq_points);

    % Estrutura para armazenar resultados de dispersão
    dispersao = struct();
    dispersao.freq_vec = freq_vec;
    dispersao.n_eff_cell = cell(1, num_freq_points);
    dispersao.beta_cell = cell(1, num_freq_points);
    dispersao.modos_por_freq = zeros(1, num_freq_points);
    
    % Opções para otimização
    options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
    
    fprintf('\n=== Iniciando cálculo de dispersão ===\n');
    fprintf('Progresso: ');

    for freq_idx = 1:num_freq_points
        
        freq = freq_vec(freq_idx);
        
        % Atualizar progresso
        if mod(freq_idx, ceil(num_freq_points/20)) == 0 || freq_idx == num_freq_points
            fprintf('%d%% ', round(freq_idx/num_freq_points*100));
        end

        % Calcular parâmetros para esta frequência
        omega = 2*pi*freq;
        k_0 = 2*pi*freq/c_0;
        
        % Definir range de busca para n_eff nesta frequência
        num_beta_points = 1000; % Pode ser reduzido para análise rápida

        n_eff_vec = linspace(n_scan_min, n_scan_max, num_beta_points);
        beta_vec = n_eff_vec * k_0;
        residuals = zeros(size(beta_vec));

        % Varredura para esta frequência
        for k = 1:num_beta_points
            residuals(k) = compute_B1(beta_vec(k), num_layers, t, layer_perms, k_0);
        end
        
        %fprintf('Refinando mínimos locais encontrados...\n');
        local_min_idxs = find(islocalmin(residuals));
        found_modes = [];
        found_betas = [];
        
        options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
        
        for idx = local_min_idxs
            beta_guess = beta_vec(idx);
            if residuals(idx) < 1.0 
                fun = @(b) compute_B1(b, num_layers, t, layer_perms, k_0);
                [beta_opt, fval] = fminsearch(fun, beta_guess, options);
                
                if fval < 1e-4 
                    n_eff_opt = beta_opt / k_0;
                    % Filtra duplicatas
                    if isempty(found_modes) || min(abs(found_modes - n_eff_opt)) > 1e-4
                        found_modes = [found_modes, n_eff_opt];
                        found_betas = [found_betas, beta_opt];
                        %hold on; plot(n_eff_opt, fval, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
                    end
                end
            end
        end

        % ORDENAR MODOS DO MAIOR PARA O MENOR n_eff
        [found_modes, sort_idx] = sort(found_modes, 'descend');
        found_betas = found_betas(sort_idx);
        
        % Armazenar resultados para esta frequência
        dispersao.n_eff_cell{freq_idx} = found_modes;
        dispersao.beta_cell{freq_idx} = found_betas;
        dispersao.modos_por_freq(freq_idx) = length(found_modes);
    end

    % ----------------------------
    % MODE-TRACKING PÓS-PROCESSO
    % ----------------------------
    % Entrada esperada:
    % dispersao.freq_vec        -> vetor de frequências (Nx1)
    % dispersao.n_eff_cell{j}   -> vetor com n_eff em cada frequência j (ordenados desc)
    % dispersao.beta_cell{j}    -> vetores de beta correspondentes
    % dispersao.modos_por_freq  -> número de modos encontrados em cada frequência
    % -------------------------------------------------------------------------
    
    freq_vec = dispersao.freq_vec;
    num_freq_points = length(freq_vec);
    
    % Determina o número máximo de modos detectados em qualquer frequência
    max_modos = max(dispersao.modos_por_freq);
    
    % Inicializa matrizes conectadas (NaN onde não houver modo)
    n_eff_conectado = NaN(num_freq_points, max_modos);
    beta_conectado  = NaN(num_freq_points, max_modos);
    
    % Marcador de candidatos já atribuídos por frequência
    assigned_idx_cell = cell(1, num_freq_points);
    for j = 1:num_freq_points
        assigned_idx_cell{j} = false(1, length(dispersao.n_eff_cell{j}));
    end
    
    % Parâmetros de tolerância (iniciais e finais)
    delta_n_abs_ini = 1e-3;
    delta_n_abs_end = 1e-5;
    
    delta_n_rel_ini = 0.02;
    delta_n_rel_end = 0.001;
    
    fprintf('Iniciando rastreamento principal de modos...\n');
    
    % =====================================================================
    % PASSO 1: Rastreamento Principal (Greedy / Semente) - Lógica Original
    % =====================================================================
    for modo = 1:max_modos
    
        % 1) Encontrar a frequência inicial onde o modo aparece
        start_freq = [];
        for j = 1:num_freq_points
            if any(~assigned_idx_cell{j})
                start_freq = j;
                break;
            end
        end
        if isempty(start_freq)
            break;  % não há mais modos
        end
    
        % 2) Escolher maior n_eff não atribuído como semente
        nlist = dispersao.n_eff_cell{start_freq};
        blist = dispersao.beta_cell{start_freq};
        free_idx = find(~assigned_idx_cell{start_freq});
    
        [~, imax_rel] = max(nlist(free_idx));
        seed_idx = free_idx(imax_rel);
    
        n_seed = nlist(seed_idx);
        b_seed = blist(seed_idx);
    
        % Atribuir semente
        n_eff_conectado(start_freq, modo) = n_seed;
        beta_conectado(start_freq, modo)  = b_seed;
        assigned_idx_cell{start_freq}(seed_idx) = true;
    
        % === Criar tolerâncias EXPO do modo ===
        freq_local = freq_vec(start_freq:end);
        w_vec = linspace(0, 1, length(freq_local));
    
        delta_n_abs_modo = delta_n_abs_ini .* (delta_n_abs_end/delta_n_abs_ini).^w_vec;
        delta_n_rel_modo = delta_n_rel_ini .* (delta_n_rel_end/delta_n_rel_ini).^w_vec;
    
        % 3) Procurar adiante
        last_n = n_seed;
        cur_freq = start_freq + 1;
    
        while cur_freq <= num_freq_points
            nlist = dispersao.n_eff_cell{cur_freq};
            blist = dispersao.beta_cell{cur_freq};
    
            if isempty(nlist)
                cur_freq = cur_freq + 1;
                continue;
            end
    
            free_idx = find(~assigned_idx_cell{cur_freq});
            if isempty(free_idx)
                cur_freq = cur_freq + 1;
                continue;
            end
    
            % Compara com last_n (valor do modo na freq anterior)
            diffs = abs(nlist(free_idx) - last_n);
            [min_diff, imin_rel] = min(diffs);
            cand_idx = free_idx(imin_rel);
            cand_n = nlist(cand_idx);
    
            % Recupera tolerância adequada
            local_idx = cur_freq - start_freq + 1;
            
            % Proteção de índice (caso local_idx exceda vetor devido a saltos)
            if local_idx > length(delta_n_abs_modo), local_idx = length(delta_n_abs_modo); end
            
            delta_n_abs = delta_n_abs_modo(local_idx);
            delta_n_rel = delta_n_rel_modo(local_idx);
    
            % Critério
            if (min_diff <= delta_n_abs) || (min_diff <= delta_n_rel * abs(last_n))
                % Atribuir
                n_eff_conectado(cur_freq, modo) = cand_n;
                beta_conectado(cur_freq, modo)  = blist(cand_idx);
                assigned_idx_cell{cur_freq}(cand_idx) = true;
    
                last_n = cand_n;
            end
            % Nota: Diferente do original, avançamos cur_freq sempre,
            % mas só atualizamos last_n se houver match.
            cur_freq = cur_freq + 1;
        end 
    end 
    
    % =====================================================================
    % PASSO 2: Varredura de Recuperação (Gap Filling)
    % =====================================================================
    fprintf('Iniciando recuperação de buracos (Gap Filling)...\n');
    
    % Tolerância "frouxa" para recuperação (evita conectar coisas absurdas)
    % 20% de variação relativa é bem grande, mas seguro para gap filling
    tol_recuperacao = 0.20; 
    
    for f = 2:num_freq_points
        
        % 1. Identificar pontos NÃO alocados nesta frequência
        idxs_sobras = find(~assigned_idx_cell{f});
        
        if isempty(idxs_sobras)
            continue; % Todos os pontos desta freq já têm dono
        end
        
        nlist_curr = dispersao.n_eff_cell{f};
        blist_curr = dispersao.beta_cell{f};
        
        % 2. Olhar para os modos conectados na frequência anterior (f-1)
        prev_modes_n = n_eff_conectado(f-1, :); % Linha com [Modo1, Modo2, ...]
        
        % Para cada ponto órfão na frequência atual
        for k = 1:length(idxs_sobras)
            orphan_idx = idxs_sobras(k);
            orphan_val = nlist_curr(orphan_idx);
            
            % Calcular distância deste órfão para TODOS os modos da freq anterior
            % Ignoramos modos que eram NaN na freq anterior
            distancias = abs(orphan_val - prev_modes_n);
            
            % Encontrar o modo da freq anterior que está mais perto deste ponto
            [min_dist, closest_mode_idx] = min(distancias);
            
            % CHECK DE SEGURANÇA:
            % A) O modo anterior existe (não é NaN)?
            if isnan(prev_modes_n(closest_mode_idx))
                continue; 
            end
            
            % B) O slot desse modo na frequência ATUAL está vazio (NaN)?
            % Se já tiver valor, o tracking principal já preencheu, então não sobrepomos.
            if ~isnan(n_eff_conectado(f, closest_mode_idx))
                continue; 
            end
            
            % C) A distância é minimamente aceitável? (Sanity check)
            if min_dist < (tol_recuperacao * abs(orphan_val))
                
                % SUCESSO: Alocar o ponto órfão ao modo encontrado
                n_eff_conectado(f, closest_mode_idx) = orphan_val;
                beta_conectado(f, closest_mode_idx)  = blist_curr(orphan_idx);
                
                % Marcar como usado (para não ser usado de novo num loop futuro se houver)
                assigned_idx_cell{f}(orphan_idx) = true;
                
                % Opcional: printar debug
                % fprintf('Recuperado ponto na freq %d para Modo %d (dist=%.4f)\n', f, closest_mode_idx, min_dist);
            end
        end
    end

    % =====================================================================
    % PASSO 3: Costura de Modos Fragmentados (Mode Stitching)
    % =====================================================================
    % Verifica se o "Modo X" que começa na freq F é na verdade
    % a continuação do "Modo Y" que terminou na freq F-1 ou F-gap.
    
    fprintf('Iniciando costura de modos fragmentados...\n');
    
    % Tolerâncias para costura
    tol_stitch_val = 0.05;  % 5% de diferença no valor do n_eff
    max_gap_size = 5;       % Máximo de frequências de buraco permitidas entre segmentos
    
    did_merge = true;
    
    % Loop "while" para garantir que se fundirmos modos, re-checamos tudo
    while did_merge
        did_merge = false;
        
        % Calcular range de existência de cada modo (Start Index e End Index)
        mode_ranges = zeros(max_modos, 2); 
        valid_modes = [];
        
        for m = 1:max_modos
            idx_validos = find(~isnan(n_eff_conectado(:, m)));
            if ~isempty(idx_validos)
                mode_ranges(m, 1) = min(idx_validos); % Start freq index
                mode_ranges(m, 2) = max(idx_validos); % End freq index
                valid_modes = [valid_modes, m];
            end
        end
        
        % Comparar todos os pares de modos válidos
        for i_idx = 1:length(valid_modes)
            mode_A = valid_modes(i_idx);
            
            for j_idx = 1:length(valid_modes)
                mode_B = valid_modes(j_idx);
                
                if mode_A == mode_B, continue; end
                
                end_A = mode_ranges(mode_A, 2);
                start_B = mode_ranges(mode_B, 1);
                
                % Checagem 1: Sequencialidade (B deve começar DEPOIS que A termina)
                % Permite um gap de até max_gap_size, mas B não pode começar antes de A acabar
                gap = start_B - end_A;
                
                if gap > 0 && gap <= max_gap_size
                    
                    % Checagem 2: Proximidade de Valor
                    val_end_A = n_eff_conectado(end_A, mode_A);
                    val_start_B = n_eff_conectado(start_B, mode_B);
                    
                    % (Opcional) Predição Linear para melhorar precisão
                    % Se A tem pelo menos 2 pontos, projeta onde estaria em start_B
                    if end_A > 1 && ~isnan(n_eff_conectado(end_A-1, mode_A))
                        slope = val_end_A - n_eff_conectado(end_A-1, mode_A); % Slope por passo
                        val_proj_A = val_end_A + slope * gap;
                    else
                        val_proj_A = val_end_A; % Sem histórico, usa valor constante
                    end
                    
                    diff_val = abs(val_start_B - val_proj_A);
                    rel_error = diff_val / abs(val_proj_A);
                    
                    if rel_error < tol_stitch_val
                        % === FUSÃO DETECTADA ===
                        fprintf('  Costurando: Modo %d (acaba freq %d) -> Modo %d (começa freq %d)\n', ...
                            mode_A, end_A, mode_B, start_B);
                        
                        % Copiar dados de B para A
                        idx_B = find(~isnan(n_eff_conectado(:, mode_B)));
                        n_eff_conectado(idx_B, mode_A) = n_eff_conectado(idx_B, mode_B);
                        beta_conectado(idx_B, mode_A)  = beta_conectado(idx_B, mode_B);
                        
                        % Apagar Modo B (preencher com NaN)
                        n_eff_conectado(:, mode_B) = NaN;
                        beta_conectado(:, mode_B)  = NaN;
                        
                        % Marcar flag para re-analisar (pois B pode ter continuidade em C)
                        did_merge = true; 
                        break; % Sai do loop interno para recalcular ranges
                    end
                end
            end
            if did_merge, break; end
        end
    end
    
    % Limpeza final: Remover colunas que ficaram vazias (NaN)
    cols_validas = any(~isnan(n_eff_conectado), 1);
    n_eff_conectado = n_eff_conectado(:, cols_validas);
    beta_conectado  = beta_conectado(:, cols_validas);
    max_modos = size(n_eff_conectado, 2); % Atualiza número real de modos


    
    
    colors = lines(max_modos);

    % ============================================================
    % PLOT ALTERNATIVO: n_eff vs EIXO X (Freq ou d/lambda_0)
    % ============================================================
    figure('Position', [100, 100, 1200, 800]);
    hold on;
    
    % Determinar qual será o eixo X
    if num_layers == 3
        % Caso com 3 camadas: Normalização d/lambda_0
        % A camada do meio é a de índice 2 (Substrato=1, Meio=2, Topo=3)
        d_meio = layer_thickness(2);
        
        % Lambda_0 para cada frequência do vetor
        lambda0_vec = c_0 ./ freq_vec;
        
        % Eixo X = d / lambda_0
        x_data_plot = d_meio ./ lambda0_vec;
        
        xlabel_str = 'Frequência Normalizada (d/\lambda_0)';
        title_str = sprintf('Dispersão Normalizada (d = %.2e m)', d_meio);
    else
        % Caso geral: Frequência em THz
        x_data_plot = freq_vec / 1e12;
        
        xlabel_str = 'Frequência [THz]';
        title_str = 'Dispersão - n_{eff} vs Frequência';
    end
    
    % Plotar cada modo conectado
    for modo_idx = 1:max_modos
        % Extrair dados para este modo
        n_eff_data = n_eff_conectado(:, modo_idx);
        
        % Remover NaNs para plotagem limpa
        valid_idx = ~isnan(n_eff_data);
        
        if any(valid_idx)
            % Seleciona apenas os pontos X e Y válidos para este modo
            x_vals = x_data_plot(valid_idx);
            y_vals = n_eff_data(valid_idx);
            
            plot(x_vals, y_vals, 'o-', ...
                'Color', colors(modo_idx, :), ...
                'MarkerSize', 4, ...
                'LineWidth', 2, ...
                'MarkerFaceColor', colors(modo_idx, :));
        end
    end
   
    % Configurar gráfico
    grid on;
    xlabel(xlabel_str, 'FontSize', 14);
    ylabel('Índice Efetivo (n_{eff})', 'FontSize', 14);
    title(title_str, 'FontSize', 14);
    set(gca, 'FontSize', 12);
   
    hold off;

    % ============================================================
    % CÁLCULO DE VELOCIDADE DE FASE E DE GRUPO
    % ============================================================
    fprintf('\n=== Calculando velocidades de fase e de grupo ===\n');
    
    % Inicializar matrizes para velocidades
    v_fase = NaN * ones(num_freq_points, max_modos);  % Velocidade de fase
    v_grupo = NaN * ones(num_freq_points, max_modos); % Velocidade de grupo
    
    % Constante da velocidade da luz no vácuo
    c = c_0; % m/s
    
    for modo_idx = 1:max_modos
        fprintf('Modo %d: calculando derivadas...\n', modo_idx);
        
        % Obter dados do modo
        n_eff_modo = n_eff_conectado(:, modo_idx);
        beta_modo = beta_conectado(:, modo_idx);
        
        % Encontrar índices válidos
        idx_validos = find(~isnan(n_eff_modo));
        
        if length(idx_validos) < 3
            fprintf('  Aviso: Poucos pontos válidos (%d) para cálculo de derivadas\n', length(idx_validos));
            continue;
        end
        
        % Para cada ponto válido
        for i = 1:length(idx_validos)
            idx = idx_validos(i);
            freq = freq_vec(idx);
            omega = 2*pi*freq;
            
            % 1. Velocidade de fase: v_phase = ω/β = c/n_eff
            v_fase(idx, modo_idx) = c / n_eff_modo(idx);
            
            % 2. Velocidade de grupo: v_group = dω/dβ
            % Precisamos calcular a derivada dβ/dω numericamente
            
            % Encontrar índices de pontos vizinhos para derivada
            idx_antes = [];
            idx_depois = [];
            
            % Procurar ponto anterior válido
            for j = idx-1:-1:max(1, idx-3)
                if ~isnan(n_eff_modo(j))
                    idx_antes = j;
                    break;
                end
            end
            
            % Procurar ponto posterior válido
            for j = idx+1:min(num_freq_points, idx+3)
                if ~isnan(n_eff_modo(j))
                    idx_depois = j;
                    break;
                end
            end
            
            % Calcular derivada usando diferentes métodos conforme disponibilidade
            if ~isempty(idx_antes) && ~isempty(idx_depois)
                % Diferença central (mais precisa)
                domega = 2*pi*(freq_vec(idx_depois) - freq_vec(idx_antes));
                dbeta = beta_modo(idx_depois) - beta_modo(idx_antes);
                v_grupo(idx, modo_idx) = domega / dbeta;
                
            elseif ~isempty(idx_antes) && isempty(idx_depois)
                % Diferença para frente (usando ponto atual e anterior)
                domega = 2*pi*(freq_vec(idx) - freq_vec(idx_antes));
                dbeta = beta_modo(idx) - beta_modo(idx_antes);
                v_grupo(idx, modo_idx) = domega / dbeta;
                
            elseif isempty(idx_antes) && ~isempty(idx_depois)
                % Diferença para trás (usando ponto atual e posterior)
                domega = 2*pi*(freq_vec(idx_depois) - freq_vec(idx));
                dbeta = beta_modo(idx_depois) - beta_modo(idx);
                v_grupo(idx, modo_idx) = domega / dbeta;
            else
                % Não há pontos suficientes para calcular derivada
                v_grupo(idx, modo_idx) = NaN;
            end
            
            % Verificar se a velocidade de grupo é fisicamente razoável
            if ~isnan(v_grupo(idx, modo_idx))
                % Velocidade de grupo não pode ser maior que c (em meios não-dispersivos)
                % mas pode ser menor que c/n_eff em meios com dispersão normal
                if v_grupo(idx, modo_idx) > 2*c || v_grupo(idx, modo_idx) < 0
                    v_grupo(idx, modo_idx) = NaN; % Valor não físico
                end
            end
        end
        
        % Suavizar velocidades de grupo (opcional)
        idx_validos_vg = find(~isnan(v_grupo(:, modo_idx)));
        if length(idx_validos_vg) > 3
            % Aplicar média móvel simples para suavizar ruído numérico
            window_size = min(3, length(idx_validos_vg));
            vg_suavizada = movmean(v_grupo(idx_validos_vg, modo_idx), window_size);
            v_grupo(idx_validos_vg, modo_idx) = vg_suavizada;
        end
    end
    

    
    % ============================================================
    % PLOT COMPARATIVO: VELOCIDADES JUNTAS POR MODO
    % ============================================================
    fprintf('\n=== Plot comparativo ===\n');
    
    % Calcular quantas figuras serão necessárias (máximo 4 modos por figura)
    modos_por_figura = 4; % 2x2
    num_figuras = ceil(max_modos / modos_por_figura);
    
    % Contador de modos já plotados
    modo_global = 0;
    
    for figura_idx = 1:num_figuras
        % Criar nova figura para cada conjunto de até 4 modos
        figure('Position', [100, 100, 1000, 800], ...
               'Name', sprintf('Velocidades - Figura %d/%d', figura_idx, num_figuras));
        
        % Para cada posição na matriz 2x2
        for posicao = 1:modos_por_figura
            modo_global = modo_global + 1;
            
            % Verificar se ainda há modos para plotar
            if modo_global > max_modos
                break;
            end
            
            % Verificar se o modo tem dados válidos
            vf_data = v_fase(:, modo_global);
            vg_data = v_grupo(:, modo_global);
            valid_idx = ~isnan(vf_data) & ~isnan(vg_data);
            
            if sum(valid_idx) < 2
                % Criar subplot vazio se não houver dados
                subplot(2, 2, posicao);
                text(0.5, 0.5, sprintf('Modo %d\nSem dados', modo_global), ...
                     'HorizontalAlignment', 'center', 'FontSize', 12);
                axis off;
                continue;
            end
            
            freqs_valid = freq_vec(valid_idx)/1e12;
            vf_valid = vf_data(valid_idx);
            vg_valid = vg_data(valid_idx);
            
            % APLICAR MÉDIA MÓVEL PARA SUAVIZAR
            janela = 5; % Tamanho da janela para média móvel
            if length(vf_valid) >= janela
                vf_suavizado = movmean(vf_valid, janela);
                vg_suavizado = movmean(vg_valid, janela);
            else
                vf_suavizado = vf_valid;
                vg_suavizado = vg_valid;
            end
            
            % Plot no mesmo gráfico com cores diferentes
            subplot(2, 2, posicao);
            hold on;
            
            plot(freqs_valid, vf_suavizado, 'b-s', ...
                'MarkerSize', 4, 'LineWidth', 1.5, ...
                'MarkerFaceColor', 'b', ...
                'DisplayName', 'v_{phase}');
            
            plot(freqs_valid, vg_suavizado, 'r-^', ...
                'MarkerSize', 4, 'LineWidth', 1.5, ...
                'MarkerFaceColor', 'r', ...
                'DisplayName', 'v_{group}');
            
            % Linha da velocidade da luz
            yline(c, '--k', 'LineWidth', 1, 'DisplayName', 'c (vácuo)');
            
            title(sprintf('Modo %d', modo_global), 'FontSize', 11, 'FontWeight', 'bold');
            xlabel('Frequência [THz]', 'FontSize', 9);
            ylabel('Velocidade [m/s]', 'FontSize', 9);
            grid on;
            
            % Formatar eixo Y
            ax = gca;
            ax.YAxis.Exponent = 6;
            
            % Legendas apenas se houver espaço
            if posicao == 1 || posicao == 3
                legend('Location', 'best', 'FontSize', 8);
            end
            
            % Ajustar limites do eixo Y para melhor visualização
            y_min = min(min(vf_suavizado), min(vg_suavizado));
            y_max = max(max(vf_suavizado), max(vg_suavizado));
            y_padding = (y_max - y_min) * 0.1;
            
            if y_padding > 0
                ylim([y_min - y_padding, y_max + y_padding]);
            end
            
            hold off;
        end
        
        fprintf('Figura %d/%d criada com modos %d-%d\n', ...
                figura_idx, num_figuras, (figura_idx-1)*modos_por_figura + 1, ...
                min(figura_idx*modos_por_figura, max_modos));
    end
    
    
    
    
    % ============================================================
    % RESUMO NUMÉRICO
    % ============================================================
    fprintf('\n=== Resumo das Velocidades ===\n');
    for modo_idx = 1:max_modos
        vf_valid = v_fase(~isnan(v_fase(:, modo_idx)), modo_idx);
        vg_valid = v_grupo(~isnan(v_grupo(:, modo_idx)), modo_idx);
        
        if ~isempty(vf_valid) && ~isempty(vg_valid)
            vf_medio = mean(vf_valid);
            vg_medio = mean(vg_valid);
            
            fprintf('Modo %d:\n', modo_idx);
            fprintf('  Velocidade de fase média: %.3e m/s (%.3f c)\n', ...
                vf_medio, vf_medio/c);
            fprintf('  Velocidade de grupo média: %.3e m/s (%.3f c)\n', ...
                vg_medio, vg_medio/c);
            fprintf('  Razão v_group/v_phase: %.3f\n', vg_medio/vf_medio);
            
            % Calcular dispersão normal/anormal
            if vg_medio < vf_medio
                fprintf('  Tipo: Dispersão normal (v_g < v_ph)\n');
            else
                fprintf('  Tipo: Dispersão anormal (v_g > v_ph)\n');
            end
        end
    end
    

    fprintf('\n=== Análise de dispersão concluída ===\n');

end


function B1 = compute_B1(beta, num_layers, t, layer_perms, k_0)
    A = zeros(num_layers,1);
    B = zeros(num_layers,1);
    % Condição de Contorno no Topo (Onda evanescente pura)
    A(num_layers) = 0;
    B(num_layers) = 1;
    
    for cont = 1:num_layers
        i = (num_layers + 1) - cont; % N, N-1 ... 1
        
        if i ~= 1
            % Calcular largura w
            % Se for a camada do topo (N), assumimos w=0 para começar na interface
            if i == num_layers
                w = 0; 
            else
                w = t(i) - t(i-1);
            end
            
            perm = layer_perms{i};
            eps_x = perm(1,1);
            eps_y = perm(3,3);
            delta = imag(perm(1,3));
            
            perm_m1 = layer_perms{i-1};
            eps_x_m1 = perm_m1(1,1);
            eps_y_m1 = perm_m1(3,3);
            delta_m1 = imag(perm_m1(1,3));
            
            gamma = sqrt(beta^2*(eps_y/eps_x) - k_0^2*((-delta^2)/eps_x + eps_y));
            theta = gamma*w;
            
            a = beta*delta - gamma*eps_x;
            b = beta*delta + gamma*eps_x;
            
            Delta_val = (-delta_m1^2 + eps_x_m1*eps_y_m1)/(-delta^2 + eps_x*eps_y);
            
            gamma_m1 = sqrt(beta^2*(eps_y_m1/eps_x_m1) - k_0^2*((-delta_m1^2)/eps_x_m1 + eps_y_m1));
            
            a_m1 = beta*delta_m1 - gamma_m1*eps_x_m1;
            b_m1 = beta*delta_m1 + gamma_m1*eps_x_m1;
            
            exp_neg = exp(-theta);
            exp_pos = exp(theta);
            
            % Propagação Matricial
            A(i-1) = A(i)*exp_neg * ((Delta_val*a - b_m1)/(a_m1 - b_m1)) + ...
                     B(i)*exp_pos  * ((Delta_val*b - b_m1)/(a_m1 - b_m1));
            
            B(i-1) = A(i)*exp_neg * ((Delta_val*a - a_m1)/(b_m1 - a_m1)) + ...
                     B(i)*exp_pos  * ((Delta_val*b - a_m1)/(b_m1 - a_m1));
        end
    end
    B1 = abs(B(1)); % O erro é o coeficiente B no substrato (que deve ser nulo)
end

function resultado_parcial_lim = resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x_j, eps_y_j, delta_j, A_j, B_j, gamma_j, t_j, lim)
% Versão robusta: trata limites r->0, casos A_j/B_j ~= 0 e evita 0/0.
% Entrada: lim pode ser limite superior/inferior (finito). Não use +/-Inf diretamente.

    % tolerâncias
    tol = 1e-12;
    exp_threshold = 700; % previne overflow em exp

    x = (lim - t_j);

    % pequenos tratamentos iniciais
    Aj = A_j; Bj = B_j;
    gj = gamma_j;
    conj_gj = conj(gj);

    r = gj + conj_gj;    % = 2*Re(gamma)
    s = gj - conj_gj;    % = 2i*Im(gamma)

    % função auxiliar: exp(div)/denom com Taylor quando denom pequeno
    function y = safeExpDiv(arg, denom)
        % arg = denom * x ou outro expoente já calculado.
        if abs(denom) < tol
            % expansão: exp(denom*x)/denom ≈ x + 0.5*denom*x^2
            y = x + 0.5*denom*(x^2);
        else
            % protege overflow/underflow
            if real(arg) > exp_threshold
                % crescimento físico -> integral diverge (retorna Inf)
                y = Inf;
            elseif real(arg) < -exp_threshold
                y = 0;
            else
                y = exp(arg)/denom;
            end
        end
    end

    % função auxiliar: exp(arg) com proteção de overflow/underflow
    function e = safeExp(arg)
        if real(arg) > exp_threshold
            e = Inf;
        elseif real(arg) < -exp_threshold
            e = 0;
        else
            e = exp(arg);
        end
    end

    % --- Termos que envolvem divisão por (gamma + conj(gamma)) ---
    % termo2: abs(A)^2 * exp((g+g*) x) / (g+g*)
    if abs(Aj) < tol
        termo2 = 0;
        termo6 = 0;
    else
        termo2 = (abs(Aj)^2) * safeExpDiv(r*x, r);
        termo6 = gj*(abs(Aj)^2) * safeExp(r*x); % usado sem divisão
    end

    % termo5: abs(B)^2 * exp(-(g+g*) x) / -(g+g*)
    if abs(Bj) < tol
        termo5 = 0;
        termo9 = 0;
    else
        termo5 = (abs(Bj)^2) * safeExpDiv((-r)*x, -r);
        termo9 = gj*(abs(Bj)^2) * safeExp((-r)*x);
    end

    % --- Termos cruzados que dividem por (g - g*) ou (g* - g) ---
    % Se s = g - g* for pequeno (i.e., g real), estes termos correspondem a 0/0.
    if abs(s) < tol
        % caso evanescente (g real) ou numericamente pequeno: os termos cruzados tendem a zero
        termo3 = 0;
        termo4 = 0;
        termo7 = 0;
        termo8 = 0;
    else
        % termo3: A*conj(B) * exp((g - g*) x)/(g - g*)
        termo3 = Aj*conj(Bj) * safeExpDiv(s*x, s);
        % termo4: B*conj(A) * exp((g - g*) x)/(conj(g) - g) = B*conj(A) * exp(s*x)/(-s)
        termo4 = Bj*conj(Aj) * safeExpDiv((-s)*x, -s);

        % termo7: gamma*conj(A)*B * exp((conj(g)-g)*x)
        termo7 = gj*(conj(Aj))*Bj * safeExp((conj_gj - gj)*x);
        % termo8: gamma*conj(B)*A * exp((g - conj(g))*x)
        termo8 = gj*(conj(Bj))*Aj * safeExp((gj - conj_gj)*x);
    end

    % Se chegamos aqui e algum termo é Inf/NaN, avisamos
    all_terms = [termo2, termo3, termo4, termo5, termo6, termo7, termo8, termo9];
    if any(isnan(all_terms)) || any(isinf(all_terms))
        warning('resultado_parcial_integral_normalizacao: termos com NaN/Inf detectados. Verifique gamma e limites.');
    end

    % agora monta os blocos conforme sua expressão original
    fator1 = 1/(2*omega*eps_0);

    termo1_prefactor = (beta * eps_y_j) / (-delta_j^2 + eps_x_j*eps_y_j);
    termo_maior1 = fator1 * ( termo1_prefactor * termo2 + termo3 + termo4 + termo5 );

    fator2 = (delta_j) / (-delta_j^2 + eps_x_j*eps_y_j);
    termo_maior2 = fator1 * ( fator2 * ( termo6 - termo7 + termo8 + termo9 ) );

    resultado_parcial_lim = termo_maior1 - termo_maior2;

    % checagens finais
    if ~isfinite(resultado_parcial_lim)
        % pode ser porque modo não é guiado (diverge) ou porque overflow
        % Retornaremos NaN e emitimos warning - o chamador deve tratar
        warning('resultado_parcial_integral_normalizacao: resultado não finito (NaN/Inf).');
    end
end


function Hy_val = calc_Hy_layer(x, A, B, gamma, t)
    % Calcula Hy para uma camada específica
    % Hy(x) = A*exp(gamma*(x-t)) + B*exp(-gamma*(x-t))
    Hy_val = A * exp(gamma*(x - t)) + B * exp(-gamma*(x - t));
end

function Ex_val = calc_Ex_layer(omega, eps_0, beta, eps_x, eps_y, delta, Hy, dHy_dx)
    % Calcula Ex para uma camada específica
    % Ex = (1/(omega*eps_0)) * (beta/(-delta^2 + eps_x*eps_y)) * (eps_y*Hy + j*delta*dHy_dx)
    fator = 1/(omega*eps_0);
    denominador = (-delta^2 + eps_x*eps_y);
    
    if abs(denominador) < 1e-12
        warning('Denominador muito pequeno em calc_Ex_layer');
        Ex_val = 0;
    else
        termo1 = (beta * eps_y) * Hy;
        termo2 = delta * dHy_dx;
        Ex_val = (fator/denominador) * (termo1 - termo2);
    end
end

function Ez_val = calc_Ez_layer(omega, eps_0, beta, eps_x, eps_y, delta, Hy, dHy_dx)
    % Calcula Ez para uma camada específica
    % Ez = (1/(omega*eps_0)) * (1/(-delta^2 + eps_x*eps_y)) * (-j*delta*eps_x*Hy + beta*eps_x*dHy_dx)
    fator = 1/(omega*eps_0);
    denominador = (-delta^2 + eps_x*eps_y);
    
    if abs(denominador) < 1e-12
        warning('Denominador muito pequeno em calc_Ez_layer');
        Ez_val = 0;
    else
        termo1 = (-1j*eps_x) * dHy_dx;
        termo2 = (1j*delta*beta) * Hy;
        Ez_val = (fator/denominador) * (termo1 + termo2);
    end
end

function dHy_dx_val = calc_dHy_dx_layer(x, A, B, gamma, t)
    % Calcula a derivada de Hy em relação a x
    % dHy/dx = gamma*A*exp(gamma*(x-t)) - gamma*B*exp(-gamma*(x-t))
    dHy_dx_val = gamma * A * exp(gamma*(x - t)) - gamma * B * exp(-gamma*(x - t));
end


function [x_points, Hy_points, Ex_points, Ez_points] = calculate_fields(num_layers, t, A, B, gamma, layer_perms, omega, eps_0, beta, num_points)
    % Calcula os campos em todo o domínio
    % num_points: número de pontos por camada para o plot
    
    % Determinar limites do domínio
    x_min = 0;
    x_max = t(end);
    
    % Criar vetor de pontos x
    total_points = num_points * num_layers;
    x_points = linspace(x_min, x_max, total_points);
    
    % Inicializar vetores de campos
    Hy_points = zeros(size(x_points));
    Ex_points = zeros(size(x_points));
    Ez_points = zeros(size(x_points));
    
    % Calcular campos para cada ponto
    for idx = 1:length(x_points)
        x = x_points(idx);
        
        % Determinar em qual camada está o ponto x
        layer_idx = find(x <= t, 1, 'first');
        if isempty(layer_idx)
            layer_idx = num_layers;
        end
        
        % Extrair parâmetros da camada
        A_layer = A(layer_idx);
        B_layer = B(layer_idx);
        gamma_layer = gamma(layer_idx);
        t_layer = t(layer_idx);
        permittivity_tensor = layer_perms{layer_idx};
        eps_x = permittivity_tensor(1,1);
        eps_y = permittivity_tensor(3,3);
        delta = imag(permittivity_tensor(1,3));
        
        % Calcular Hy e sua derivada
        Hy_val = calc_Hy_layer(x, A_layer, B_layer, gamma_layer, t_layer);
        dHy_dx_val = calc_dHy_dx_layer(x, A_layer, B_layer, gamma_layer, t_layer);
        
        % Calcular Ex e Ez
        Ex_val = calc_Ex_layer(omega, eps_0, beta, eps_x, eps_y, delta, Hy_val, dHy_dx_val);
        Ez_val = calc_Ez_layer(omega, eps_0, beta, eps_x, eps_y, delta, Hy_val, dHy_dx_val);
        
        % Armazenar resultados
        Hy_points(idx) = Hy_val;
        Ex_points(idx) = Ex_val;
        Ez_points(idx) = Ez_val;
    end
end


function plot_fields(modo_idx, n_eff, x_points, Hy_points, Ex_points, Ez_points, t, num_layers, size_order_layers)
    % Plotar os campos com interfaces das camadas
    
    % Converter para unidades apropriadas para plot
    x_plot = x_points / size_order_layers;
    t_plot = t / size_order_layers;
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot para Hy (Campo Magnético)
    subplot(3,1,1);
    plot(x_plot, real(Hy_points), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Re(Hy)');
    hold on;
    plot(x_plot, imag(Hy_points), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Im(Hy)');
    plot(x_plot, abs(Hy_points), 'g--', 'LineWidth', 1, 'DisplayName', '|Hy|');
    
    % Adicionar linhas verticais para interfaces
    for i = 1:length(t_plot)
        xline(t_plot(i), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end

    ylabel('H_y [A/m]');
    title('H_y');
    legend('Location', 'best');
    grid on;
    
    % Subplot para Ex (Campo Elétrico)
    subplot(3,1,2);
    plot(x_plot, real(Ex_points), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Re(Ex)');
    hold on;
    plot(x_plot, imag(Ex_points), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Im(Ex)');
    plot(x_plot, abs(Ex_points), 'g--', 'LineWidth', 1, 'DisplayName', '|Ex|');
    
    % Adicionar linhas verticais para interfaces
    for i = 1:length(t_plot)
        xline(t_plot(i), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end
    
    ylabel('E_x [V/m]');
    title('E_x');
    legend('Location', 'best');
    grid on;
    
    % Subplot para Ez (Campo Elétrico)
    subplot(3,1,3);
    plot(x_plot, real(Ez_points), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Re(Ez)');
    hold on;
    plot(x_plot, imag(Ez_points), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Im(Ez)');
    plot(x_plot, abs(Ez_points), 'g--', 'LineWidth', 1, 'DisplayName', '|Ez|');
    
    % Adicionar linhas verticais para interfaces
    for i = 1:length(t_plot)
        xline(t_plot(i), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end
    
    xlabel(['Posição x [', get_unit_label(size_order_layers), ']']);
    ylabel('E_z [V/m]');
    title('E_z');
    legend('Location', 'best');
    grid on;
    
    % Adicionar título geral
    sgtitle(sprintf('Campos do Modo %d (n_{eff} = %.4f)', modo_idx, n_eff));
end

function unit_label = get_unit_label(size_order)
    % Retorna o label da unidade baseado na ordem de grandeza
    if size_order == 1e-9
        unit_label = 'nm';
    elseif size_order == 1e-6
        unit_label = 'μm';
    elseif size_order == 1e-3
        unit_label = 'mm';
    else
        unit_label = 'm';
    end
end

function [A, B, gamma, fator_normalizacao, resultado_integral] = processar_modos_encontrados(beta, num_layers, t, layer_perms, eps_0, omega, k_0)
% PROCESSAR_MODOS_ENCONTRADOS Processa cada modo encontrado calculando coeficientes e normalizando
%
% Entradas:
%   found_modes     - Vetor de índices efetivos encontrados
%   found_betas     - Vetor de constantes de propagação correspondentes
%   num_layers      - Número de camadas da estrutura
%   t               - Vetor de posições das interfaces
%   layer_perms     - Celula com tensores de permissividade de cada camada
%   eps_0           - Permissividade do vácuo
%   omega           - Frequência angular
%   k_0             - Número de onda no vácuo
%
% Saídas:
%   resultados_modos - Estrutura com resultados para cada modo
%   A_array         - Matriz com coeficientes A para cada modo (colunas)
%   B_array         - Matriz com coeficientes B para cada modo (colunas)
%   gamma_array     - Matriz com parâmetros gamma para cada modo (colunas)

    fprintf('\n--- Processando cada modo encontrado ---\n');
        
        % Inicializar vetores para este modo
        A = zeros(num_layers, 1);
        B = zeros(num_layers, 1);
        gamma = zeros(num_layers, 1);
        w = zeros(num_layers, 1);
        theta = zeros(num_layers, 1);
        a = zeros(num_layers, 1);
        b = zeros(num_layers, 1);
        
        % Condição de contorno no topo
        A(num_layers) = 0;
        B(num_layers) = 1;
        
        resultado_integral = 0;
        
        % Propagação para obter coeficientes A e B
        for cont = 1:num_layers
            i = (num_layers + 1) - cont;
            
            if i ~= 1
                w(i) = t(i) - t(i - 1);
                
                % Parâmetros da camada i
                permittivity_tensor = layer_perms{i};
                eps_x = permittivity_tensor(1,1);
                eps_y = permittivity_tensor(3,3);
                delta = imag(permittivity_tensor(1,3));
                
                % Parâmetros da camada i-1
                permittivity_tensor_minus1 = layer_perms{i-1};
                eps_x_minus1 = permittivity_tensor_minus1(1,1);
                eps_y_minus1 = permittivity_tensor_minus1(3,3);
                delta_minus1 = imag(permittivity_tensor_minus1(1,3));
                
                % Cálculo dos parâmetros de propagação
                gamma(i) = sqrt(beta^2 * (eps_y/eps_x) - k_0^2 * ((-delta^2)/eps_x + eps_y));
                
  
                % Correção física: a parte real de gamma deve ser positiva para representar decaimento
                if real(gamma(i)) < 0
                    gamma(i) = -gamma(i);
                end
                
                theta(i) = w(i) * gamma(i);
                a(i) = beta * delta - gamma(i) * eps_x;
                b(i) = beta * delta + gamma(i) * eps_x;
                
                Delta_val = (-delta_minus1^2 + eps_x_minus1 * eps_y_minus1) / ...
                           (-delta^2 + eps_x * eps_y);
                
                gamma_minus1 = sqrt(beta^2 * (eps_y_minus1/eps_x_minus1) - ...
                                   k_0^2 * ((-delta_minus1^2)/eps_x_minus1 + eps_y_minus1));
                a_minus1 = beta * delta_minus1 - gamma_minus1 * eps_x_minus1;
                b_minus1 = beta * delta_minus1 + gamma_minus1 * eps_x_minus1;
                
                % Propagação matricial
                A(i-1) = A(i) * exp(-theta(i)) * ((Delta_val * a(i) - b_minus1)/(a_minus1 - b_minus1)) + ...
                        B(i) * exp(theta(i)) * ((Delta_val * b(i) - b_minus1)/(a_minus1 - b_minus1));
                
                B(i-1) = A(i) * exp(-theta(i)) * ((Delta_val * a(i) - a_minus1)/(b_minus1 - a_minus1)) + ...
                        B(i) * exp(theta(i)) * ((Delta_val * b(i) - a_minus1)/(b_minus1 - a_minus1));
                
              
            else
                % Última camada (substrato)
                w(i) = t(i);
                
                permittivity_tensor = layer_perms{i};
                eps_x = permittivity_tensor(1,1);
                eps_y = permittivity_tensor(3,3);
                delta = imag(permittivity_tensor(1,3));
                
                gamma(i) = sqrt(beta^2 * (eps_y/eps_x) - k_0^2 * ((-delta^2)/eps_x + eps_y));
                
            end
        end
        
        resultado_integral = 0;
    
        for cont = 1:num_layers
            i = (num_layers + 1) - cont;
            
            permittivity_tensor = layer_perms{i};
            eps_x = permittivity_tensor(1,1);
            eps_y = permittivity_tensor(3,3);
            delta = imag(permittivity_tensor(1,3));
        
            if i ~= 1
                resultado_parcial = resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A(i), B(i), gamma(i), t(i), t(i)) - ...
                                   resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A(i), B(i), gamma(i), t(i), t(i-1));
            else
                resultado_parcial = resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A(i), B(i), gamma(i), t(i), t(i)) - ...
                                   resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A(i), B(i), gamma(i), t(i), 0);
            end
            
            resultado_integral = resultado_integral + resultado_parcial;
        end
        
        fprintf('Integral ANTES da normalização: %.12g W/m\n', real(resultado_integral));
        
        % 3. TERCEIRO: Calcular fator de normalização
        if abs(resultado_integral) > 1e-20
            fator_normalizacao = 1/sqrt(resultado_integral);
        else
            warning('Integral de normalização muito pequena. Verifique se o modo é guiado.');
            fator_normalizacao = 1;
        end
        
        fprintf('Fator de normalização: %.12g\n', fator_normalizacao);
        
        % 4. QUARTO: Aplicar normalização a TODOS os coeficientes A e B
        A_normalized = A * fator_normalizacao;
        B_normalized = B * fator_normalizacao;
        
        % 5. QUINTO: VERIFICAR a normalização (deve resultar em 1 W/m)
        resultado_integral_normalizacao = 0;
        
        for cont = 1:num_layers
            i = (num_layers + 1) - cont;
            
            permittivity_tensor = layer_perms{i};
            eps_x = permittivity_tensor(1,1);
            eps_y = permittivity_tensor(3,3);
            delta = imag(permittivity_tensor(1,3));
        
            if i == 1
                resultado_parcial = resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A_normalized(i), B_normalized(i), gamma(i), t(i), t(i)) - ...
                                   resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A_normalized(i), B_normalized(i), gamma(i), t(i), 0);
            else
                resultado_parcial = resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A_normalized(i), B_normalized(i), gamma(i), t(i), t(i)) - ...
                                   resultado_parcial_integral_normalizacao(eps_0, omega, beta, eps_x, eps_y, delta, A_normalized(i), B_normalized(i), gamma(i), t(i), t(i-1));
            end
            
            resultado_integral_normalizacao = resultado_integral_normalizacao + resultado_parcial;
        end
        
        fprintf('Integral DEPOIS da normalização: %.12g W/m\n', real(resultado_integral_normalizacao));
        fprintf('Diferença de 1 W/m: %.2e W/m\n', abs(real(resultado_integral_normalizacao) - 1));
        
        teste = abs(real(resultado_integral_normalizacao) - 1) < 1e-6;
        if teste
            fprintf('✓ Normalização bem-sucedida! Integral = 1 W/m (precisão 1e-6)\n');
        else
            fprintf('⚠ Atenção: Normalização não atingiu precisão desejada\n');
        end
        
        % Atribuir os valores normalizados de volta às variáveis principais
        A = A_normalized;
        B = B_normalized;

end

function B1_complex = compute_B1_muller(beta, num_layers, t, layer_perms, k_0)
    A = zeros(num_layers,1);
    B = zeros(num_layers,1);
    
    % Condição de Contorno no Topo
    A(num_layers) = 0;
    B(num_layers) = 1;
    
    for cont = 1:num_layers
        i = (num_layers + 1) - cont; 
        
        if i ~= 1
            if i == num_layers, w = 0; else, w = t(i) - t(i-1); end
            
            perm = layer_perms{i};
            eps_x = perm(1,1); eps_y = perm(3,3); delta = imag(perm(1,3));
            perm_m1 = layer_perms{i-1};
            eps_x_m1 = perm_m1(1,1); eps_y_m1 = perm_m1(3,3); delta_m1 = imag(perm_m1(1,3));
            
            % --- CORREÇÃO FUNDAMENTAL AQUI ---
            gamma = sqrt(beta^2*(eps_y/eps_x) - k_0^2*((-delta^2)/eps_x + eps_y));
            if real(gamma) < 0, gamma = -gamma; end % Força decaimento físico
            % ---------------------------------
            
            theta = gamma*w;
            a = beta*delta - gamma*eps_x;
            b = beta*delta + gamma*eps_x;
            
            Delta_val = (-delta_m1^2 + eps_x_m1*eps_y_m1)/(-delta^2 + eps_x*eps_y);
            
            % --- CORREÇÃO TAMBÉM NA CAMADA ANTERIOR ---
            gamma_m1 = sqrt(beta^2*(eps_y_m1/eps_x_m1) - k_0^2*((-delta_m1^2)/eps_x_m1 + eps_y_m1));
            if real(gamma_m1) < 0, gamma_m1 = -gamma_m1; end % Força decaimento físico
            % ------------------------------------------
            
            a_m1 = beta*delta_m1 - gamma_m1*eps_x_m1;
            b_m1 = beta*delta_m1 + gamma_m1*eps_x_m1;
            
            exp_neg = exp(-theta);
            exp_pos = exp(theta);
            
            A(i-1) = A(i)*exp_neg * ((Delta_val*a - b_m1)/(a_m1 - b_m1)) + ...
                     B(i)*exp_pos  * ((Delta_val*b - b_m1)/(a_m1 - b_m1));
            
            B(i-1) = A(i)*exp_neg * ((Delta_val*a - a_m1)/(b_m1 - a_m1)) + ...
                     B(i)*exp_pos  * ((Delta_val*b - a_m1)/(b_m1 - a_m1));
        end
    end
    B1_complex = B(1); 
end

% ============================================================
% ALGORITMO DE MULLER (DO SEU ARQUIVO ENVIADO)
% ============================================================
function f_val = muller(f, x0, x1, x2)
    % Function implements Muller's method
    iter_max = 100; % max number of steps in Muller method
    f_tol = 1e-9;   % Tolerância ajustada
    x_tol = 1e-9;
    
    y0 = f(x0);
    y1 = f(x1);
    y2 = f(x2);
    
    iter = 0;
    while(iter <= iter_max)
        iter = iter + 1;
        
        % Cálculo dos coeficientes da parábola
        q = (x1 - x2) / (x1 - x0); % dummy var to avoid singular matrix if possible? 
        % Usando a forma explícita padrão do Muller (mais robusta):
        h1 = x1 - x0;
        h2 = x2 - x1;
        d1 = (y1 - y0) / h1;
        d2 = (y2 - y1) / h2;
        
        a = (d2 - d1) / (h2 + h1);
        b = a*h2 + d2;
        c = y2;
        
        if (a~=0)
            D = sqrt(b*b - 4*a*c);
            % Escolhe o sinal para maximizar denominador (menor passo)
            if abs(b - D) < abs(b + D)
                E = b + D;
            else
                E = b - D;
            end
            dx = -2*c / E;
        elseif (b~=0)
            dx = -c/b;
        else
            dx = 0; % Falha ou convergencia plana
        end
        
        x3 = x2 + dx;
        
        % Atualização dos pontos
        x0 = x1; y0 = y1;
        x1 = x2; y1 = y2;
        x2 = x3; y2 = f(x2);
        
        % Critério de parada
        if (abs(dx) < x_tol || abs(y2) < f_tol)
            break;
        end
    end
    
    if (abs(y2) < 1e-4) % Verifica se convergiu mesmo
        f_val = x2;
    else
        f_val = 0;
    end
end

function [X_mesh, Y_mesh, Field_2D] = expand_fields_to_xy(x_points, field_1D, width_y, num_points_y)
    % Expande o campo 1D (x) para uma matriz 2D (x,y) assumindo invariância em y.
    %
    % Entradas:
    %   x_points: vetor de posições x (vertical/espessura)
    %   field_1D: vetor de campo complexo calculado (ex: Hy_points)
    %   width_y: largura arbitrária para visualização em y (ex: 5 microns)
    %   num_points_y: resolução na direção y
    
    % Criar vetor Y (centralizado em 0)
    y_vec = linspace(-width_y/2, width_y/2, num_points_y);
    
    % Criar a malha (Meshgrid)
    % Note que X varia nas linhas e Y nas colunas para facilitar a visualização matricial
    [Y_mesh, X_mesh] = meshgrid(y_vec, x_points);
    
    % Replicar o campo 1D para todas as colunas de Y
    % repmat repete o vetor coluna field_1D, 1 vez na vertical, N vezes na horizontal
    Field_2D = repmat(field_1D(:), 1, num_points_y);
end

function plot_2D_heatmaps(modo_idx, n_eff, X_mesh, Y_mesh, Hy_2D, Ex_2D, Ez_2D, t, size_order_layers)
    % Plota os heatmaps dos campos no plano transversal XY
    
    % Conversão de unidades para labels
    unit_str = get_unit_label(size_order_layers);
    
    % Normalizar coordenadas para plotagem (tirar da ordem de grandeza, ex: micrômetros)
    X_plot = X_mesh / size_order_layers;
    Y_plot = Y_mesh / size_order_layers;
    t_plot = t / size_order_layers;
    
    % Criar figura
    figure('Position', [150, 150, 1400, 500], 'Name', sprintf('Modo %d - Heatmaps 2D (XY)', modo_idx));
    
    % --- Plot Hy ---
    subplot(1, 3, 1);
    % Usamos abs() para mostrar a magnitude/intensidade do campo
    pcolor(Y_plot, X_plot, abs(Hy_2D)); 
    shading interp; % Suaviza as cores (remove linhas de grade)
    colormap('jet');
    colorbar;
    hold on;
    % Desenhar interfaces das camadas
    for k = 1:length(t_plot)
        yline(t_plot(k), 'w--', 'LineWidth', 1.5); % Linha branca tracejada
    end
    title('|H_y| (Magnitude)');
    xlabel(['Largura y [' unit_str ']']);
    ylabel(['Espessura x [' unit_str ']']);
    axis tight;
    
    % --- Plot Ex ---
    subplot(1, 3, 2);
    pcolor(Y_plot, X_plot, abs(Ex_2D)); 
    shading interp;
    colormap('jet');
    colorbar;
    hold on;
    for k = 1:length(t_plot)
        yline(t_plot(k), 'w--', 'LineWidth', 1.5);
    end
    title('|E_x| (Magnitude)');
    xlabel(['Largura y [' unit_str ']']);
    % ylabel não precisa repetir
    axis tight;
    
    % --- Plot Ez ---
    subplot(1, 3, 3);
    pcolor(Y_plot, X_plot, abs(Ez_2D)); 
    shading interp;
    colormap('jet');
    colorbar;
    hold on;
    for k = 1:length(t_plot)
        yline(t_plot(k), 'w--', 'LineWidth', 1.5);
    end
    title('|E_z| (Magnitude)');
    xlabel(['Largura y [' unit_str ']']);
    axis tight;
    
    sgtitle(sprintf('Distribuição Transversal de Campos (Modo %d, n_{eff}=%.4f)', modo_idx, real(n_eff)));
end