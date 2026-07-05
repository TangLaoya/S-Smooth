function test_smoothing_diagnostics
% 诊断脚本：验证Z和S平滑效果
% 对比Test_beforeFit.s26p与output.s26p，分析是否存在跳变

close all;
clear;

n = 26;
snp = ['s' num2str(n) 'p'];
inputFile  = ['Test_beforeFit.' snp];
outputFile = ['output.' snp];

fprintf('========== 平滑效果诊断 ==========\n\n');

% 1. 读取输入和输出文件
[freq_in, S_in, params_in] = readSnP_simple(inputFile);
[freq_out, S_out, params_out] = readSnP_simple(outputFile);

[Nf, Nport, ~] = size(S_in);
fprintf('输入文件：频点数=%d，端口数=%d\n', Nf, Nport);

% 2. 转换为Z
Z0 = 50;
I_mat = eye(Nport);
Z_in = zeros(Nf, Nport, Nport);
Z_out = zeros(Nf, Nport, Nport);

for k = 1:Nf
    Sk = squeeze(S_in(k,:,:));
    Z_in(k,:,:) = S_to_Z_stable(Sk, Z0, I_mat);
    
    Sk = squeeze(S_out(k,:,:));
    Z_out(k,:,:) = S_to_Z_stable(Sk, Z0, I_mat);
end

% 3. 检测平滑段起点
M_opt = detectWidestSmoothSegment(freq_in, S_in, Nport);
fprintf('平滑频段起点: %d (频率 %.4e Hz)\n', M_opt, freq_in(M_opt));

% 4. 对每个端口对分析Z曲线平滑性
fprintf('\n========== Z曲线平滑性分析 ==========\n');

% 提取低频段（M_opt前面的部分）
low_freq_idx = 1:min(M_opt+10, Nf);  % 包括平滑段起点之后的一些点
freq_low = freq_in(low_freq_idx);
logf_low = log10(freq_low);

% 统计信息
jump_info = [];
smooth_improvement = [];

for i = 1:Nport
    for j = 1:Nport
        Z_curve_in = squeeze(Z_in(low_freq_idx, i, j));
        Z_curve_out = squeeze(Z_out(low_freq_idx, i, j));
        
        % 计算幅值曲线的差分
        amp_in = abs(Z_curve_in);
        amp_out = abs(Z_curve_out);
        
        amp_in(amp_in < 1e-30) = 1e-30;
        amp_out(amp_out < 1e-30) = 1e-30;
        
        dB_in = 20*log10(amp_in);
        dB_out = 20*log10(amp_out);
        
        % 计算对数坐标下的差分（衡量平滑性）
        diff_in_dB = abs(diff(dB_in));
        diff_out_dB = abs(diff(dB_out));
        
        % 最大跳变
        max_jump_in = max(diff_in_dB);
        max_jump_out = max(diff_out_dB);
        
        % 平均跳变
        mean_jump_in = mean(diff_in_dB);
        mean_jump_out = mean(diff_out_dB);
        
        % 判断是否仍有跳变（超过2dB）
        has_jump_out = any(diff_out_dB > 2.0);
        
        if has_jump_out || (i == 1 && j == 1)  % 只输出前几个或有跳变的
            fprintf('Z(%d,%d): 输入[%.3f, %.3f]dB → 输出[%.3f, %.3f]dB (max, mean jump)\n', ...
                i, j, max_jump_in, mean_jump_in, max_jump_out, mean_jump_out);
        end
        
        % 记录跳变信息
        if has_jump_out
            jump_info = [jump_info; i, j, max_jump_out, mean_jump_out];
        end
        
        if max_jump_in > 0.1
            improvement = (max_jump_in - max_jump_out) / max_jump_in * 100;
            smooth_improvement = [smooth_improvement; improvement];
        end
    end
end

if isempty(jump_info)
    fprintf('\n✓ 所有Z曲线在低频段都达到平滑（最大跳变 < 2.0 dB）\n');
else
    fprintf('\n⚠ 检测到 %d 个Z曲线仍存在跳变（>2.0 dB）\n', size(jump_info, 1));
    fprintf('  建议进一步优化外推和插值算法\n');
end

% 平均改善
if ~isempty(smooth_improvement)
    fprintf('平均改善率: %.2f%%\n\n', mean(smooth_improvement(smooth_improvement > 0)));
end

% 5. 对数坐标下的线性拟合检验
fprintf('========== 对数坐标线性拟合检验 ==========\n');

% 计算低频段在对数坐标下的残差
linearity_error = [];

for i = 1:min(Nport, 3)  % 只检查前3个端口对
    for j = 1:min(Nport, 3)
        Z_curve_out = squeeze(Z_out(low_freq_idx, i, j));
        amp_out = abs(Z_curve_out);
        amp_out(amp_out < 1e-30) = 1e-30;
        dB_out = 20*log10(amp_out);
        
        % 线性拟合
        p = polyfit(logf_low, dB_out, 1);
        fit_vals = polyval(p, logf_low);
        residuals = dB_out - fit_vals;
        rms_error = sqrt(mean(residuals.^2));
        
        linearity_error = [linearity_error; rms_error];
        
        fprintf('Z(%d,%d): 线性拟合RMS误差 = %.4f dB\n', i, j, rms_error);
    end
end

if ~isempty(linearity_error)
    fprintf('平均线性度误差: %.4f dB\n\n', mean(linearity_error));
end

% 6. 绘图
fprintf('========== 生成诊断图表 ==========\n');

% 选择前3个有跳变的端口对或前3个端口对
plot_idx = [];
if ~isempty(jump_info)
    plot_idx = jump_info(1:min(3, size(jump_info,1)), 1:2);
else
    for i = 1:min(2, Nport)
        for j = 1:min(2, Nport)
            plot_idx = [plot_idx; i, j];
            if size(plot_idx, 1) >= 4, break; end
        end
        if size(plot_idx, 1) >= 4, break; end
    end
end

if ~isempty(plot_idx)
    fig = figure('NumberTitle', 'off', 'Name', 'Z-Parameter平滑性诊断');
    n_plots = size(plot_idx, 1);
    
    for p = 1:n_plots
        i = plot_idx(p, 1);
        j = plot_idx(p, 2);
        
        Z_curve_in = squeeze(Z_in(low_freq_idx, i, j));
        Z_curve_out = squeeze(Z_out(low_freq_idx, i, j));
        
        amp_in = abs(Z_curve_in);
        amp_out = abs(Z_curve_out);
        
        amp_in(amp_in < 1e-30) = 1e-30;
        amp_out(amp_out < 1e-30) = 1e-30;
        
        dB_in = 20*log10(amp_in);
        dB_out = 20*log10(amp_out);
        
        subplot(n_plots, 1, p);
        loglog(freq_low, amp_in, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Input');
        hold on;
        loglog(freq_low, amp_out, 's-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Output');
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('|Z| (Ω)');
        title(sprintf('Z(%d,%d) - Low Frequency Segment', i, j));
        legend('Location', 'best');
        hold off;
    end
    
    savefig(fig, 'smoothing_diagnostics.fig');
    fprintf('诊断图已保存为: smoothing_diagnostics.fig\n');
end

% 7. 输出建议
fprintf('\n========== 诊断建议 ==========\n');

if isempty(jump_info)
    fprintf('✓ 平滑效果良好，Z曲线在低频段保持连续平滑\n');
    fprintf('✓ 可以使用输出文件进行后续处理\n');
else
    fprintf('⚠ 平滑效果仍需改进：\n');
    fprintf('  1. 考虑调整外推权重或约束阈值\n');
    fprintf('  2. 检查M_opt（平滑段起点）的探测是否准确\n');
    fprintf('  3. 可能需要增加插值点的密度\n');
    fprintf('  4. 考虑使用样条插值而非线性插值\n');
end

fprintf('\n========== 诊断完成 ==========\n');

end

% ======================= 函数复制（与主程序一致）=======================

function Z = S_to_Z_stable(S, Z0, I_mat)
    denom = I_mat - S;
    cond_num = cond(denom);
    if cond_num > 1e12
        denom = denom + 1e-10 * I_mat;
    end
    numer = I_mat + S;
    Z = Z0 * (numer / denom);
end

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
