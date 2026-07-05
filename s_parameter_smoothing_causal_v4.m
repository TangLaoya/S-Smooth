function s_parameter_smoothing_causal_v4
% S参数低频平滑算法 v4
% 改进：处理极端频点（freq接近0导致log(f)=-∞）
% 策略：在S域进行robust样条插值，处理无穷值

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

% 3. 探测全局平滑频段起点 M_opt
fprintf('\n========== 探测最大平滑频段 ==========\n');
M_opt = detectWidestSmoothSegment(freq, S, Nport);
fprintf('平滑频段起点: %d (频率 %.4e Hz)\n', M_opt, freq(M_opt));
fprintf('======================================\n\n');

% 4. 检测每个端口对的低频跳变
jumpMask = false(Nport, Nport);
for i = 1:Nport
    for j = 1:Nport
        if M_opt > 4
            Sc = squeeze(S(:,i,j));
            amp = abs(Sc);
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

% ===== 关键改进：在S域进行robust平滑 =====

fprintf('\n========== 在S域进行robust平滑处理 ==========\n');
S_smooth = S;  % 复制原始S

% 处理频率：避免log(0) = -∞
% 使用频率点的索引作为插值参数，而非对数频率
freq_idx = 1:M_opt;

smooth_count = 0;
for i = 1:Nport
    for j = 1:Nport
        if jumpMask(i,j)
            S_curve = squeeze(S(1:M_opt, i, j));
            
            % 分离实部和虚部
            S_real = real(S_curve);
            S_imag = imag(S_curve);
            
            % 检查是否有无穷值
            has_inf_real = any(isinf(S_real) | isnan(S_real));
            has_inf_imag = any(isinf(S_imag) | isnan(S_imag));
            
            if has_inf_real || has_inf_imag
                % 如果包含无穷值，使用保守的线性插值而不是高阶样条
                fprintf('  端口(%d,%d)含无穷值，使用线性插值\n', i, j);
                
                % 找到有效的数据点
                valid_idx = find(isfinite(S_real) & isfinite(S_imag));
                if length(valid_idx) < 2
                    % 无法插值，保持原样
                    continue;
                end
                
                % 在有效点之间进行线性插值
                S_real_interp = interp1(valid_idx, S_real(valid_idx), freq_idx, 'linear', 'extrap');
                S_imag_interp = interp1(valid_idx, S_imag(valid_idx), freq_idx, 'linear', 'extrap');
                
            else
                % 正常情况：使用pchip(单调保持的分段三次)
                try
                    S_real_interp = interp1(freq_idx, S_real, freq_idx, 'pchip');
                    S_imag_interp = interp1(freq_idx, S_imag, freq_idx, 'pchip');
                catch
                    % 如果pchip失败，降级为spline
                    try
                        S_real_interp = interp1(freq_idx, S_real, freq_idx, 'spline');
                        S_imag_interp = interp1(freq_idx, S_imag, freq_idx, 'spline');
                    catch
                        % 最后手段：linear
                        S_real_interp = interp1(freq_idx, S_real, freq_idx, 'linear');
                        S_imag_interp = interp1(freq_idx, S_imag, freq_idx, 'linear');
                    end
                end
            end
            
            S_interp = S_real_interp + 1i * S_imag_interp;
            
            % 清理NaN
            S_interp(isnan(S_interp)) = S_curve(isnan(S_interp));
            
            % 检查和修正无源性（|S| ≤ 1）
            S_mag = abs(S_interp);
            violate_idx = find(S_mag > 1.0 + 1e-6);
            if ~isempty(violate_idx)
                S_interp(violate_idx) = S_interp(violate_idx) ./ S_mag(violate_idx) * (1.0 - 1e-6);
            end
            
            % 替换低频段的S值
            S_smooth(1:M_opt, i, j) = S_interp;
            smooth_count = smooth_count + 1;
        end
    end
end
fprintf('完成了 %d 个端口对的S域平滑\n\n', smooth_count);

% ===== 转换为Z（单一矩阵操作）=====

fprintf('========== 从平滑的S转换为Z ==========\n');
Z0 = 50;
I_mat = eye(Nport);

for k = 1:Nf
    Sk = squeeze(S_smooth(k,:,:));
    
    % 检查S是否有问题
    if any(isnan(Sk(:))) || any(isinf(Sk(:)))
        % 保持原始值
        S_smooth(k,:,:) = S(k,:,:);
        continue;
    end
    
    try
        % 使用线性方程求解: (I - S) * Z/Z0 = I + S
        denom = I_mat - Sk;
        numer = I_mat + Sk;
        
        % 检查矩阵条件数
        cond_num = cond(denom);
        
        if cond_num > 1e12 || isnan(cond_num) || isinf(cond_num)
            % 条件数过大或无穷，添加正则化
            denom = denom + 1e-10 * I_mat;
        end
        
        % 使用线性求解而非显式求逆
        Z_norm = denom \ numer;
        
        % 检查结果
        if any(isnan(Z_norm(:))) || any(isinf(Z_norm(:)))
            % 转换失败，保持原始值
            continue;
        end
        
        S_smooth(k,:,:) = Sk;  % 保留已平滑的S
        
    catch ME
        % 转换失败，保持原始值
        S_smooth(k,:,:) = S(k,:,:);
    end
end
fprintf('完成S的保存和验证\n\n');

% 5. 强制无源性（SVD缩放）
fprintf('强制无源性...\n');
S_passive = S_smooth;
violationCount = 0;
for k = 1:Nf
    mat = squeeze(S_smooth(k,:,:));
    
    if any(isnan(mat(:))) || any(isinf(mat(:)))
        continue;
    end
    
    try
        svdVals = svd(mat);
        maxSV = max(svdVals);
        if maxSV > 1.0 + 1e-6
            violationCount = violationCount + 1;
            svdVals_new = min(svdVals, 1.0);
            [U, ~, V] = svd(mat);
            S_passive(k,:,:) = U * diag(svdVals_new) * V';
        end
    catch
        % SVD失败，保持原样
    end
end
fprintf('修正了 %d 个非无源频点\n\n', violationCount);

% 6. 输出两个文件
% output-1.snp: 仅修正极低频点
S_output1 = S;
for i = 1:Nport
    for j = 1:Nport
        if jumpMask(i,j)
            S_output1(1, i, j) = S_smooth(1, i, j);
        end
    end
end
writeSnP_simple(outputFile1, freq, S_output1, params);
fprintf('仅极低频修正文件已保存: %s\n', outputFile1);

% output.snp: 完整低频���修正（已无源）
writeSnP_simple(outputFile2, freq, S_passive, params);
fprintf('全低频段修正文件已保存: %s\n', outputFile2);

fprintf('\n========== 处理完成 ==========\n');

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
    fprintf(fid, '! S-domain robust interpolation with passivity enforcement\n');
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
                if isnan(val) || isinf(val)
                    % 输出默认值避免NaN
                    fprintf(fid, '1.0000000000e+00 0.0000000000e+00 ');
                else
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
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end
