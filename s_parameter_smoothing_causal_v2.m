function s_parameter_smoothing_causal
% Z域实部/虚部物理外推替换 + 无源性强制（SVD缩放）
% 改进版：Z→S转换严格约束 + 数值稳定性优化 + 对数频率插值

n = 26;
snp = ['s' num2str(n) 'p'];
inputFile  = ['Test_beforeFit.' snp];
outputFile1 = ['output-1.' snp];
outputFile2 = ['output.' snp];

% 1. 读取
[freq, S, params] = readSnP_simple(inputFile);
[Nf, Nport, ~] = size(S);
fprintf('读取成功：端口数=%d，实际频点数=%d\n', Nport, Nf);

% 2. 预处理
[freq, S] = preprocessFreq(freq, S);
S(isnan(S)) = 1e-15 + 1i*1e-15;
S(isinf(S)) = 1e10 + 1i*1e10;

% 3. 转换为Z（数值稳定版）
Z0 = 50;
I_mat = eye(Nport);
Z = zeros(Nf, Nport, Nport);
for k = 1:Nf
    Sk = squeeze(S(k,:,:));
    Z(k,:,:) = S_to_Z_stable(Sk, Z0, I_mat);
end

% 4. 探测全局平滑频段起点 M_opt
fprintf('\n========== 探测最大平滑频段 ==========\n');
M_opt = detectWidestSmoothSegment(freq, S, Nport);
fprintf('平滑频段起点: %d (频率 %.4e Hz)\n', M_opt, freq(M_opt));
fprintf('======================================\n\n');

% 5. 检测每个端口对的低频跳变
jumpMask = false(Nport, Nport);
for i = 1:Nport
    for j = 1:Nport
        if M_opt > 4
            Zc = squeeze(Z(:,i,j));
            amp = abs(Zc);
            amp(amp < 1e-30) = 1e-30;
            dB = 20*log10(amp);
            logf = log10(freq);
            x_low = logf(1:M_opt-1);
            y_low = dB(1:M_opt-1);
            valid = isfinite(x_low) & isfinite(y_low);
            if sum(valid) < 5
                continue;
            end
            p = polyfit(x_low(valid), y_low(valid), 1);
            resid = std(y_low(valid) - polyval(p, x_low(valid)));
            if resid > 1.5
                dev = abs(y_low(valid) - polyval(p, x_low(valid)));
                count = 0;
                for k = 1:length(dev)
                    if dev(k) > 1.5
                        count = count + 1;
                        if count >= 3
                            jumpMask(i,j) = true;
                            break;
                        end
                    else
                        count = 0;
                    end
                end
            end
        end
    end
end
fprintf('检测到 %d 个低频跳变端口对\n', sum(jumpMask(:)));

% ===== 新的外推与插值逻辑 =====

% 6. 改进的外推：同时处理整个矩阵，保证S基本不变
fprintf('\n========== 改进的Z矩阵外推（S约束）==========\n');
Z_corrected = Z;
S_low_corrected = zeros(M_opt-1, Nport, Nport);

for k_low = 1:(M_opt-1)
    f_curr = freq(k_low);
    f_ref = freq(M_opt);
    
    for i = 1:Nport
        for j = 1:Nport
            if jumpMask(i,j)
                % 获取参考频点（平滑段起始）的Z和S
                Z_ref = squeeze(Z(M_opt, i, j));
                S_original = S(k_low, i, j);
                
                % 第一阶段：基于物理规律的Z外推
                [Z_extrap_init, extrap_type] = extrapolate_Z_physical(f_curr, f_ref, Z_ref, Z0);
                
                % 第二阶段：检查S变化，若过大则调整R
                S_extrap_init = Z_to_S_element(Z_extrap_init, Z0);
                S_delta = abs(S_extrap_init - S_original);
                
                if S_delta > 0.05 && abs(S_original) > 0.01  % S变化超过阈值
                    % 调用约束优化：最小化S变化，同时保持X的物理特性
                    Z_extrap = constrained_Z_adjustment(Z_extrap_init, S_original, Z0, ...
                                                         f_curr, f_ref, extrap_type);
                else
                    Z_extrap = Z_extrap_init;
                end
                
                S_low_corrected(k_low, i, j) = Z_to_S_element(Z_extrap, Z0);
                Z_corrected(k_low, i, j) = Z_extrap;
            else
                % 非跳变元素保持原样
                S_low_corrected(k_low, i, j) = S(k_low, i, j);
                Z_corrected(k_low, i, j) = Z(k_low, i, j);
            end
        end
    end
end

fprintf('完成外推和S约束调整\n\n');

% 7. 改进的插值：在对数频率空间进行线性插值
fprintf('========== 对数频率空间线性插值 ==========\n');
logf = log10(freq);

for i = 1:Nport
    for j = 1:Nport
        if jumpMask(i,j)
            % 被插值的两个端点：最低频点和平滑段起始点
            Z_low = squeeze(Z_corrected(1, i, j));
            Z_ref = squeeze(Z_corrected(M_opt, i, j));
            
            % 对数频率
            logf_low = logf(1);
            logf_ref = logf(M_opt);
            
            % 在对数频率空间中对Z的实部和虚部分别进行线性插值
            for k = 1:M_opt
                % 插值权重（在对数频率空间）
                w = (logf(k) - logf_low) / (logf_ref - logf_low);
                w = max(0, min(1, w));  % 限制在[0,1]
                
                % 线性插值
                R_interp = (1-w) * real(Z_low) + w * real(Z_ref);
                X_interp = (1-w) * imag(Z_low) + w * imag(Z_ref);
                
                Z_corrected(k, i, j) = R_interp + 1i * X_interp;
            end
        end
    end
end

fprintf('完成对数频率空间线性插值\n\n');

% 8. 修正后的Z转回S（使用稳定转换）
S_corrected = zeros(Nf, Nport, Nport);
for k = 1:Nf
    Zk = squeeze(Z_corrected(k,:,:));
    S_corrected(k,:,:) = Z_to_S_matrix_stable(Zk, Z0, I_mat);
end

% 9. 强制无源性（SVD缩放）
fprintf('强制无源性...\n');
S_passive = S_corrected;
violationCount = 0;
for k = 1:Nf
    mat = squeeze(S_corrected(k,:,:));
    svdVals = svd(mat);
    maxSV = max(svdVals);
    if maxSV > 1.0 + 1e-6
        violationCount = violationCount + 1;
        svdVals_new = min(svdVals, 1.0);
        [U, ~, V] = svd(mat);
        S_passive(k,:,:) = U * diag(svdVals_new) * V';
    end
end
fprintf('修正了 %d 个非无源频点\n', violationCount);

% 10. 输出两个文件
% output-1.snp: 仅修正极低频点（第一个频点）
S_output1 = S;
for i = 1:Nport
    for j = 1:Nport
        if jumpMask(i,j)
            S_output1(1, i, j) = S_low_corrected(1, i, j);
        end
    end
end
writeSnP_simple(outputFile1, freq, S_output1, params);
fprintf('仅极低频修正文件已保存: %s\n', outputFile1);

% output.snp: 完整低频段修正（已无源）
writeSnP_simple(outputFile2, freq, S_passive, params);
fprintf('全低频段修正文件已保存: %s\n', outputFile2);

end

% ======================= 新增核心函数 =======================

function Z = S_to_Z_stable(S, Z0, I_mat)
    % 数值稳定的S→Z单矩阵变换
    % Z = Z0 * (I + S) * inv(I - S)
    denom = I_mat - S;
    
    % 计算条件数检测奇异性
    cond_num = cond(denom);
    if cond_num > 1e12
        % 使用正则化
        denom = denom + 1e-10 * I_mat;
    end
    
    numer = I_mat + S;
    
    % 使用更稳定的求解方式而不是显式求逆
    Z = Z0 * (numer / denom);
end

function S = Z_to_S_matrix_stable(Z, Z0, I_mat)
    % 数值稳定的Z→S矩阵变换
    % S = (Z/Z0 - I) / (Z/Z0 + I)
    z_norm = Z / Z0;
    
    numer = z_norm - I_mat;
    denom = z_norm + I_mat;
    
    % 计算条件数检测奇异性
    cond_num = cond(denom);
    if cond_num > 1e12
        % 使用正则化
        denom = denom + 1e-10 * I_mat;
    end
    
    S = numer / denom;
end

function S_elem = Z_to_S_element(Z_elem, Z0)
    % 单个Z元素→S元素的转换
    % S = (Z/Z0 - 1) / (Z/Z0 + 1)
    z_norm = Z_elem / Z0;
    denom = z_norm + 1;
    
    if abs(denom) < 1e-15
        S_elem = 1.0 + 1e-15i;
    else
        S_elem = (z_norm - 1) / denom;
    end
end

function [Z_extrap, type] = extrapolate_Z_physical(f_curr, f_ref, Z_ref, Z0)
    % 基于物理规律的Z外推
    % 返回外推的Z值和外推类型
    
    R_ref = real(Z_ref);
    X_ref = imag(Z_ref);
    
    % 判断元素类型
    if abs(X_ref) < 1e-9 * (abs(R_ref) + 1e-12)
        type = 'resistive';
        R_curr = R_ref;
        X_curr = 0;
    elseif X_ref < 0
        type = 'capacitive';
        R_curr = R_ref;
        % 容性：X ∝ 1/f，低频时X更大（负向）
        X_curr = X_ref * (f_ref / f_curr);
    else
        type = 'inductive';
        R_curr = R_ref;
        % 感性：X ∝ f，低频时X更小
        X_curr = X_ref * (f_curr / f_ref);
    end
    
    Z_extrap = R_curr + 1i * X_curr;
end

function Z_adj = constrained_Z_adjustment(Z_init, S_target, Z0, f_curr, f_ref, extrap_type)
    % 约束优化：调整Z以使S基本保持不变
    % 关键约束：最小化 |Z_to_S_element(Z_adj, Z0) - S_target|
    
    R_init = real(Z_init);
    X_init = imag(Z_init);
    
    % 使用分析性方法：对于高阻抗，S对Z实部敏感度低
    % 我们采用迭代校正：主要调整R使得S匹配
    
    S_init = Z_to_S_element(Z_init, Z0);
    S_error = S_init - S_target;
    
    % 计算Z对S的导数（近似）
    dZ_small = 0.1 * (1 + abs(R_init));  % 扰动幅度
    Z_perturb_R = (R_init + dZ_small) + 1i*X_init;
    S_perturb_R = Z_to_S_element(Z_perturb_R, Z0);
    dS_dR = (S_perturb_R - S_init) / dZ_small;
    
    % 计算所需的R调整（一阶泰勒）
    if abs(dS_dR) > 1e-15
        delta_R = -S_error / dS_dR * 0.5;  % 保守因子0.5防止过度校正
        R_adj = R_init + delta_R;
    else
        R_adj = R_init;
    end
    
    % 虚部保持外推值不变（物理意义强）
    Z_adj = R_adj + 1i * X_init;
    
    % 再次检查S误差
    S_adj = Z_to_S_element(Z_adj, Z0);
    S_error_new = abs(S_adj - S_target);
    
    if S_error_new > abs(S_error)  % 校正变差，回退
        Z_adj = Z_init;
    end
end

% ======================= 辅助函数 =======================

function M_opt = detectWidestSmoothSegment(freq, S, Nport)
    Nf = length(freq);
    smooth_mask_all = false(Nf, Nport, Nport);
    for i = 1:Nport
        for j = 1:Nport
            curve = squeeze(S(:,i,j));
            smooth_mask_all(:,i,j) = extractSmoothMask(curve, Nf);
        end
    end
    smooth_ratio = sum(smooth_mask_all, [2,3]) / (Nport*Nport);
    threshold_ratio = 0.5;
    potential_idx = find(smooth_ratio >= threshold_ratio);
    if ~isempty(potential_idx)
        [longest_start, longest_end] = findLongestContinuousSegment(smooth_ratio >= threshold_ratio);
        if (longest_end - longest_start + 1 < 0.2*Nf) || (longest_end < 0.7*Nf)
            longest_start = round(Nf * 0.4);
            longest_end = Nf;
        else
            if longest_end - longest_start + 1 > 0.9*Nf
                longest_start = round(Nf * 0.3);
            end
        end
    else
        longest_start = round(Nf * 0.4);
        longest_end = Nf;
    end
    if longest_end - longest_start < 20
        longest_start = max(1, Nf - 30);
        longest_end = Nf;
    end
    M_opt = longest_start;
end

function mask = extractSmoothMask(curve, Nf)
    amp = abs(curve);
    amp(amp < 1e-30) = 1e-30;
    log_amp = 20*log10(amp);
    diff_log = abs(diff(log_amp));
    med_val = median(diff_log);
    mad_val = median(abs(diff_log - med_val));
    if mad_val < 1e-15
        mad_val = std(diff_log);
    end
    threshold = med_val + 8 * mad_val + 1e-12;
    mask = [true; diff_log < threshold];
    for k = 2:Nf-1
        if mask(k) && mask(k-1) && mask(k+1)
            mask(k) = true;
        else
            mask(k) = false;
        end
    end
    if sum(mask) < 5
        threshold = med_val + 12 * mad_val + 1e-12;
        mask = [true; diff_log < threshold];
        for k = 2:Nf-1
            if mask(k) && mask(k-1) && mask(k+1)
                mask(k) = true;
            else
                mask(k) = false;
            end
        end
    end
    mask = mask(:);
end

function [start_idx, end_idx] = findLongestContinuousSegment(mask)
    mask = mask(:);
    if ~any(mask)
        start_idx = 1; end_idx = 1;
        return;
    end
    diff_mask = diff([0; mask; 0]);
    starts = find(diff_mask == 1);
    ends = find(diff_mask == -1) - 1;
    if isempty(starts)
        start_idx = 1; end_idx = 1;
        return;
    end
    lengths = ends - starts + 1;
    [~, idx] = max(lengths);
    start_idx = starts(idx);
    end_idx = ends(idx);
end

% ======================= I/O 辅助函数 =======================
function [freq, S, params] = readSnP_simple(filename)
    fid = fopen(filename, 'r');
    assert(fid ~= -1, '无法打开文件: %s', filename);
    wholeFile = fread(fid, '*char')';
    fclose(fid);
    hashLine = regexp(wholeFile, '#[^\n\r]*', 'match', 'once');
    if isempty(hashLine), error('未找到 # 选项行'); end
    tokens = strsplit(strtrim(hashLine));
    freq_unit = upper(tokens{2});
    param_type = tokens{3};
    format_type = upper(tokens{4});
    R0 = 50;
    for k = 5:numel(tokens)-1
        if strcmpi(tokens{k}, 'R')
            R0 = str2double(tokens{k+1}); break;
        end
    end
    Nport_declared = str2double(regexp(wholeFile, '\[Number of Ports\]\s*(\d+)', 'tokens', 'once'));
    startIdx = strfind(wholeFile, '[Network Data]');
    if isempty(startIdx), error('缺少 [Network Data]'); end
    dataStr = wholeFile(startIdx+14:end);
    endIdx = strfind(dataStr, '[End of Data]');
    if ~isempty(endIdx), dataStr(endIdx:end) = []; end
    dataStr = regexprep(dataStr, '!.*?\n', ' ');
    allNums = sscanf(dataStr, '%f');
    totalNums = length(allNums);
    if isempty(Nport_declared)
        Nport = sqrt((totalNums - 1)/2);
        if mod(Nport,1)~=0, error('无法推断端口数'); end
    else
        Nport = Nport_declared;
    end
    inferredNf = floor(totalNums / (1 + 2*Nport^2));
    if inferredNf < 1, error('数据量不足'); end
    Nfreq = inferredNf;
    validNums = Nfreq * (1 + 2*Nport^2);
    allNums = allNums(1:validNums);
    scale = getScale(freq_unit);
    freq = zeros(Nfreq,1); S = zeros(Nfreq,Nport,Nport);
    idx = 1;
    for k = 1:Nfreq
        freq(k) = allNums(idx) * scale;
        idx = idx+1;
        for i = 1:Nport
            for j = 1:Nport
                a = allNums(idx); b = allNums(idx+1);
                switch format_type
                    case 'RI', S(k,i,j) = a + 1i*b;
                    case 'MA', S(k,i,j) = a.*exp(1i*b*pi/180);
                    case 'DB', S(k,i,j) = 10^(a/20).*exp(1i*b*pi/180);
                end
                idx = idx+2;
            end
        end
    end
    params = struct('R0',R0,'param_type',param_type,'format_type',format_type,'freq_unit',freq_unit,'scale',scale);
end

function scl = getScale(unit)
    switch upper(unit)
        case 'HZ', scl=1;
        case 'KHZ', scl=1e3;
        case 'MHZ', scl=1e6;
        case 'GHZ', scl=1e9;
        otherwise, error('未知频率单位');
    end
end

function [freq, S] = preprocessFreq(freq, S)
    [freqSorted, order] = sort(freq);
    S = S(order,:,:);
    [~, ia] = unique(freqSorted, 'stable');
    dupIdx = setdiff(1:numel(freqSorted), ia);
    if ~isempty(dupIdx)
        fprintf('检测到 %d 个重复频率点，已取平均。\n', numel(dupIdx));
        [freq, ~, ic] = unique(freqSorted);
        Nf = numel(freq);
        S_new = zeros(Nf, size(S,2), size(S,3));
        for k = 1:Nf
            mask = (ic == k);
            S_new(k,:,:) = mean(S(mask,:,:), 1);
        end
        S = S_new;
    else
        freq = freqSorted;
    end
end

function writeSnP_simple(filename, freq, S, params)
    fid = fopen(filename, 'w');
    Nport = size(S,2); Nf = size(S,1);
    fprintf(fid, '# %s %s %s R %.1f\n', params.freq_unit, params.param_type, params.format_type, params.R0);
    fprintf(fid, '! Z-domain R/X physical extrapolation with S-parameter constraint and passivity enforcement\n');
    fprintf(fid, '! Ports: %d, Frequencies: %d\n', Nport, Nf);
    outScale = 1;
    switch upper(params.freq_unit)
        case 'KHZ', outScale=1e-3;
        case 'MHZ', outScale=1e-6;
        case 'GHZ', outScale=1e-9;
    end
    for k = 1:Nf
        fprintf(fid, '%.8e ', freq(k)*outScale);
        for i = 1:Nport
            for j = 1:Nport
                val = S(k,i,j);
                switch params.format_type
                    case 'RI'
                        fprintf(fid, '%.10e %.10e ', real(val), imag(val));
                    case 'MA'
                        fprintf(fid, '%.10e %.10e ', abs(val), angle(val)*180/pi);
                    case 'DB'
                        db = -1e3;
                        if abs(val)>0, db=20*log10(abs(val)); end
                        fprintf(fid, '%.10e %.10e ', db, angle(val)*180/pi);
                end
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end
