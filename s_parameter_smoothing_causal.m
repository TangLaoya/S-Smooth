function s_parameter_smoothing_causal
% Z域实部/虚部物理外推替换 + 无源性强制（SVD缩放）

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

% 3. 转换为Z（加正则化）
Z0 = 50;
I_mat = eye(Nport); % 定义单位矩阵
Z = zeros(Nf, Nport, Nport);
for k = 1:Nf
    Sk = squeeze(S(k,:,:));
    Zk = Z0 * (I_mat + Sk) / (I_mat - Sk + 1e-12*I_mat);
    Z(k,:,:) = Zk;
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

% 6. 对跳变端口对进行实部/虚部物理外推替换
Z_corrected = Z;
for i = 1:Nport
    for j = 1:Nport
        if jumpMask(i,j)
            Zc = squeeze(Z(:,i,j));
            
            % 锚点：平滑段起始点
            f_ref = freq(M_opt);
            Z_ref = Zc(M_opt);
            R_dc = real(Z_ref);
            X_ref = imag(Z_ref);
            
            Z_new_low = zeros(M_opt-1, 1);
            for k = 1:M_opt-1
                f_curr = freq(k);
                
                % 实部外推：低频趋于直流电阻常数
                R_curr = R_dc;
                
                % 虚部外推：基于物理规律
                if abs(X_ref) < 1e-9 * (abs(R_dc) + 1e-12)
                    % 纯电阻
                    X_curr = 0;
                elseif X_ref < 0
                    % 容性：X = -1/(2*pi*f*C)，与频率成反比，低频趋于 -Inf
                    X_curr = X_ref * (f_ref / f_curr);
                else
                    % 感性：X = 2*pi*f*L，与频率成正比，低频趋于 0
                    X_curr = X_ref * (f_curr / f_ref);
                end
                
                Z_new_low(k) = R_curr + 1i * X_curr;
            end
            
            % 替换低频段
            Z_corrected(1:M_opt-1, i, j) = Z_new_low;
        end
    end
end

% 7. 修正后的Z转回S（加正则化，修复原代码I未定义的bug）
S_corrected = zeros(Nf, Nport, Nport);
for k = 1:Nf
    Zk = squeeze(Z_corrected(k,:,:));
    % 修复：使用 I_mat 替代未定义的 I，并统一正则化项
    Sk = (Zk/Z0 - I_mat) / (Zk/Z0 + I_mat + 1e-12*I_mat);
    S_corrected(k,:,:) = Sk;
end

% 8. 强制无源性（SVD缩放）
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

% 9. 输出两个文件
% output-1.snp: 仅修正极低频点（第一个频点）
S_output1 = S;
for i = 1:Nport
    for j = 1:Nport
        if jumpMask(i,j)
            Zk_dc = Z_corrected(1,i,j);
            % 修复转S计算
            S_output1(1,i,j) = (Zk_dc/Z0 - 1) / (Zk_dc/Z0 + 1);
        end
    end
end
writeSnP_simple(outputFile1, freq, S_output1, params);
fprintf('仅极低频修正文件已保存: %s\n', outputFile1);

% output.snp: 完整低频段修正（已无源）
writeSnP_simple(outputFile2, freq, S_passive, params);
fprintf('全低频段修正文件已保存: %s\n', outputFile2);

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
    fprintf(fid, '! Z-domain R/X physical extrapolation with passivity enforcement\n');
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