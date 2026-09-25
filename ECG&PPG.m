%% =========================================================
% ECG & PPG Analysis
%
% Channel_1 = ECG
% Channel_2 = PPG
%
% 功能：
% 1. 一次選擇一個或多個 TXT
% 2. 自動由 Time 計算取樣頻率
% 3. ECG 5-15 Hz band-pass filter
% 4. 自動偵測 R wave
% 5. 計算 ECG Heart Rate
% 6. PPG 0.5-8 Hz band-pass filter
% 7. 自動偵測 PPG systolic peak
% 8. 計算 PPG Pulse Rate
% 9. 計算 PTT
% 10. 偵測重搏切跡 dicrotic notch 候選點
% 11. 偵測 notch 後的 secondary peak
% 12. 顯示各種偵測結果
% 13. 輸出總整理表
%
% 需要：
% Signal Processing Toolbox
%
% 使用函數：
% butter
% filtfilt
% findpeaks
%% =========================================================

clear;
clc;
close all;


%% =========================================================
% 1. 選擇資料
%% =========================================================

[fileNames, filePath] = uigetfile( ...
    '*.txt', ...
    '選擇 ECG & PPG TXT 資料', ...
    'MultiSelect', 'on');

if isequal(fileNames, 0)
    disp('沒有選擇檔案');
    return;
end

% 如果只選一個檔案，轉成 cell
if ischar(fileNames)
    fileNames = {fileNames};
end

nFiles = length(fileNames);


%% =========================================================
% 2. 建立結果儲存空間
%% =========================================================

Result_File = strings(nFiles,1);

Result_fs = zeros(nFiles,1);
Result_Duration = zeros(nFiles,1);

Result_HR = zeros(nFiles,1);
Result_PulseRate = zeros(nFiles,1);

Result_PTT_mean = zeros(nFiles,1);
Result_PTT_std = zeros(nFiles,1);
Result_PTT_median = zeros(nFiles,1);

Result_Rnumber = zeros(nFiles,1);
Result_PTTnumber = zeros(nFiles,1);

% 重搏切跡
Result_NotchNumber = zeros(nFiles,1);
Result_NotchRate = zeros(nFiles,1);


%% =========================================================
% 3. 分析每一個檔案
%% =========================================================

for k = 1:nFiles

    fileName = fileNames{k};
    fullName = fullfile(filePath, fileName);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('分析檔案：%s\n', fileName);
    fprintf('============================================\n');


    %% -----------------------------------------------------
    % 3.1 讀取 TXT
    %% -----------------------------------------------------

    fid = fopen(fullName, 'r');

    if fid == -1
        error('無法開啟檔案：%s', fileName);
    end

    % 跳過第一列標題
    fgetl(fid);

    % 三欄：
    % Time / Channel_1 / Channel_2
    data = textscan(fid, '%f %f %f');

    fclose(fid);

    t   = data{1};
    ECG = data{2};
    PPG = data{3};


    %% -----------------------------------------------------
    % 3.2 基本資料檢查
    %% -----------------------------------------------------

    if isempty(t)
        warning('%s 沒有讀到資料', fileName);
        continue;
    end

    dt = median(diff(t));

    fs = 1 / dt;

    duration = t(end) - t(1);

    fprintf('取樣頻率 = %.1f Hz\n', fs);
    fprintf('記錄時間 = %.2f s\n', duration);
    fprintf('樣本數   = %d\n', length(t));


    %% =====================================================
    % 4. 去除 DC offset
    %% =====================================================

    ECG0 = ECG - mean(ECG);
    PPG0 = PPG - mean(PPG);


    %% =====================================================
    % 5. ECG Band-pass Filter
    %
    % 5-15 Hz 用於突顯 QRS
    %% =====================================================

    [bECG, aECG] = butter( ...
        3, ...
        [5 15] / (fs/2), ...
        'bandpass');

    ECG_f = filtfilt( ...
        bECG, ...
        aECG, ...
        ECG0);


    %% =====================================================
    % 6. R-wave Detection
    %% =====================================================

    % 最小 R-R 間距 0.4 秒
    % 約等於最大心率 150 bpm
    minRR = round(0.40 * fs);

    % 自動 threshold
    R_prominence = 1.0 * std(ECG_f);

    [R_amp, R_loc] = findpeaks( ...
        ECG_f, ...
        'MinPeakDistance', minRR, ...
        'MinPeakProminence', R_prominence);

    R_time = t(R_loc);


    %% =====================================================
    % 7. ECG Heart Rate
    %% =====================================================

    RR = diff(R_time);

    % 合理 R-R interval：
    % 0.4-1.5 秒
    validRR = ...
        RR > 0.40 & ...
        RR < 1.50;

    RR_valid = RR(validRR);

    if isempty(RR_valid)

        HR = NaN;

    else

        HR = 60 / mean(RR_valid);

    end

    fprintf('\nECG：\n');
    fprintf('R wave 數量 = %d\n', length(R_loc));
    fprintf('平均心率 = %.1f bpm\n', HR);


    %% =====================================================
    % 8. PPG Band-pass Filter
    %
    % 0.5-8 Hz 保留主要脈波形態
    %% =====================================================

    [bPPG, aPPG] = butter( ...
        3, ...
        [0.5 8] / (fs/2), ...
        'bandpass');

    PPG_f = filtfilt( ...
        bPPG, ...
        aPPG, ...
        PPG0);


    %% =====================================================
    % 9. PPG Systolic Peak Detection
    %% =====================================================

    minPP = round(0.40 * fs);

    PPG_prominence = ...
        0.30 * std(PPG_f);

    [PPG_amp, PPG_loc] = findpeaks( ...
        PPG_f, ...
        'MinPeakDistance', minPP, ...
        'MinPeakProminence', PPG_prominence);

    PPG_time = t(PPG_loc);


    %% =====================================================
    % 10. PPG Pulse Rate
    %% =====================================================

    PP_interval = diff(PPG_time);

    validPP = ...
        PP_interval > 0.40 & ...
        PP_interval < 1.50;

    PP_valid = PP_interval(validPP);

    if isempty(PP_valid)

        PulseRate = NaN;

    else

        PulseRate = 60 / mean(PP_valid);

    end

    fprintf('\nPPG：\n');
    fprintf('PPG 主峰數量 = %d\n', length(PPG_loc));
    fprintf('平均脈率 = %.1f bpm\n', PulseRate);


    %% =====================================================
    % 11. PTT Calculation
    %
    % PTT =
    % PPG systolic peak time - ECG R-wave time
    %
    % 尋找每個 R wave 後
    % 80-600 ms 內第一個 PPG peak
    %% =====================================================

    PTT = [];

    paired_R = [];
    paired_PPG = [];

    for i = 1:length(R_time)

        candidate = find( ...
            PPG_time > R_time(i) + 0.080 & ...
            PPG_time < R_time(i) + 0.600);

        if ~isempty(candidate)

            j = candidate(1);

            currentPTT = ...
                PPG_time(j) - R_time(i);

            PTT(end+1,1) = ...
                currentPTT;

            paired_R(end+1,1) = ...
                R_time(i);

            paired_PPG(end+1,1) = ...
                PPG_time(j);

        end

    end


    %% =====================================================
    % 12. PTT Outlier Removal
    %
    % 使用 Median Absolute Deviation
    %% =====================================================

    if length(PTT) >= 3

        medPTT = median(PTT);

        MAD_PTT = ...
            median(abs(PTT - medPTT));

        if MAD_PTT > 0

            sigma_est = ...
                1.4826 * MAD_PTT;

            validPTT = ...
                abs(PTT - medPTT) <= ...
                3 * sigma_est;

            PTT_clean = ...
                PTT(validPTT);

        else

            PTT_clean = PTT;

        end

    else

        PTT_clean = PTT;

    end


    %% =====================================================
    % 13. PTT Statistics
    %% =====================================================

    if isempty(PTT_clean)

        PTT_mean = NaN;
        PTT_std = NaN;
        PTT_median = NaN;

    else

        % 秒 -> ms
        PTT_mean = ...
            mean(PTT_clean) * 1000;

        PTT_std = ...
            std(PTT_clean) * 1000;

        PTT_median = ...
            median(PTT_clean) * 1000;

    end

    fprintf('\nPTT：\n');

    fprintf( ...
        '有效 PTT 數量 = %d\n', ...
        length(PTT_clean));

    fprintf( ...
        '平均 PTT = %.1f ms\n', ...
        PTT_mean);

    fprintf( ...
        'PTT 標準差 = %.1f ms\n', ...
        PTT_std);

    fprintf( ...
        'PTT 中位數 = %.1f ms\n', ...
        PTT_median);


    %% =====================================================
    % 14. Dicrotic Notch Detection
    %
    % 原則：
    %
    % systolic peak
    %       |
    %       v
    % 主峰 -> notch valley -> secondary peak
    %
    % 搜尋範圍：
    % 主峰後 80-400 ms
    %
    % 注意：
    % 這裡屬於候選點偵測，
    % 最後仍需人工確認波形
    %% =====================================================

    notch_loc = [];
    secondary_loc = [];

    notch_time = [];
    secondary_time = [];

    % 主峰後 80 ms 開始搜尋
    notchStartDelay = ...
        round(0.08 * fs);

    % 最晚搜尋到 400 ms
    notchEndDelay = ...
        round(0.40 * fs);


    for i = 1:length(PPG_loc)

        systolicPeak = ...
            PPG_loc(i);

        startIndex = ...
            systolicPeak + notchStartDelay;

        endIndex = ...
            systolicPeak + notchEndDelay;


        % 如果有下一個 pulse peak
        % 不讓搜尋範圍跨到下一個心動週期
        if i < length(PPG_loc)

            nextPeakLimit = ...
                PPG_loc(i+1) - round(0.05*fs);

            endIndex = ...
                min(endIndex, nextPeakLimit);

        end


        % 避免超過訊號長度
        if endIndex > length(PPG_f)

            endIndex = length(PPG_f);

        end


        if startIndex >= endIndex

            continue;

        end


        % ---------------------------------
        % 找 valley
        % ---------------------------------

        segment = ...
            PPG_f(startIndex:endIndex);

        notchProminence = ...
            0.03 * std(PPG_f);

        [~, valleyLoc] = findpeaks( ...
            -segment, ...
            'MinPeakProminence', ...
            notchProminence);


        if isempty(valleyLoc)

            continue;

        end


        % ---------------------------------
        % 逐個 valley 測試
        % 後方是否存在 secondary peak
        % ---------------------------------

        notchFound = false;

        for v = 1:length(valleyLoc)

            currentNotch = ...
                startIndex + valleyLoc(v) - 1;


            % notch 後至少 20 ms
            secondaryStart = ...
                currentNotch + round(0.02*fs);

            % notch 後最多搜尋 250 ms
            secondaryEnd = ...
                currentNotch + round(0.25*fs);


            % 不跨到下一個 systolic peak
            if i < length(PPG_loc)

                secondaryEnd = ...
                    min(secondaryEnd, ...
                    PPG_loc(i+1) - 1);

            end


            secondaryEnd = ...
                min(secondaryEnd, ...
                length(PPG_f));


            if secondaryStart >= secondaryEnd

                continue;

            end


            secondarySegment = ...
                PPG_f(secondaryStart:secondaryEnd);


            secondaryProminence = ...
                0.03 * std(PPG_f);


            [~, secLoc] = findpeaks( ...
                secondarySegment, ...
                'MinPeakProminence', ...
                secondaryProminence);


            if isempty(secLoc)

                continue;

            end


            % 找到 notch 後的第一個 secondary peak
            currentSecondary = ...
                secondaryStart + secLoc(1) - 1;


            notch_loc(end+1,1) = ...
                currentNotch;

            secondary_loc(end+1,1) = ...
                currentSecondary;

            notchFound = true;

            break;

        end

    end


    % 轉成時間
    if ~isempty(notch_loc)

        notch_time = ...
            t(notch_loc);

        secondary_time = ...
            t(secondary_loc);

    end


    %% =====================================================
    % 15. Dicrotic Notch Statistics
    %% =====================================================

    NotchNumber = ...
        length(notch_loc);

    if isempty(PPG_loc)

        NotchRate = NaN;

    else

        NotchRate = ...
            NotchNumber / length(PPG_loc) * 100;

    end


    fprintf('\n重搏切跡：\n');

    fprintf( ...
        '偵測到重搏切跡候選點 = %d\n', ...
        NotchNumber);

    fprintf( ...
        '重搏切跡出現比例 = %.1f %%\n', ...
        NotchRate);


    %% =====================================================
    % 16. 畫圖
    %% =====================================================

    figure( ...
        'Name', fileName, ...
        'NumberTitle', 'off');


    % -----------------------------------------------------
    % Raw ECG
    % -----------------------------------------------------

    subplot(4,1,1);

    plot(t, ECG);

    ylabel('ECG');

    title( ...
        ['Raw ECG - ', fileName], ...
        'Interpreter', 'none');

    grid on;


    % -----------------------------------------------------
    % Raw PPG
    % -----------------------------------------------------

    subplot(4,1,2);

    plot(t, PPG);

    ylabel('PPG');

    title('Raw PPG');

    grid on;


    % -----------------------------------------------------
    % ECG + R wave
    % -----------------------------------------------------

    subplot(4,1,3);

    plot(t, ECG_f);

    hold on;

    plot( ...
        R_time, ...
        ECG_f(R_loc), ...
        'ro', ...
        'MarkerSize', 5);

    ylabel('Filtered ECG');

    title( ...
        sprintf( ...
        'R-wave Detection   HR = %.1f bpm', ...
        HR));

    legend( ...
        'Filtered ECG', ...
        'R wave');

    grid on;


    % -----------------------------------------------------
    % PPG + systolic peak + notch + secondary peak
    % -----------------------------------------------------

    subplot(4,1,4);

    plot(t, PPG_f);

    hold on;


    % PPG systolic peak
    plot( ...
        PPG_time, ...
        PPG_f(PPG_loc), ...
        'ro', ...
        'MarkerSize', 5);


    % Dicrotic notch
    if ~isempty(notch_loc)

        plot( ...
            notch_time, ...
            PPG_f(notch_loc), ...
            'bs', ...
            'MarkerSize', 5);

    end


    % Secondary peak
    if ~isempty(secondary_loc)

        plot( ...
            secondary_time, ...
            PPG_f(secondary_loc), ...
            'g^', ...
            'MarkerSize', 5);

    end


    xlabel('Time (s)');
    ylabel('Filtered PPG');

    title( ...
        sprintf( ...
        'Pulse = %.1f bpm   PTT = %.1f ms   Notch = %.1f %%', ...
        PulseRate, ...
        PTT_mean, ...
        NotchRate));


    legend( ...
        'Filtered PPG', ...
        'Systolic Peak', ...
        'Dicrotic Notch', ...
        'Secondary Peak');


    grid on;


    %% =====================================================
    % 17. 儲存結果
    %% =====================================================

    Result_File(k) = ...
        string(fileName);

    Result_fs(k) = ...
        fs;

    Result_Duration(k) = ...
        duration;

    Result_HR(k) = ...
        HR;

    Result_PulseRate(k) = ...
        PulseRate;

    Result_PTT_mean(k) = ...
        PTT_mean;

    Result_PTT_std(k) = ...
        PTT_std;

    Result_PTT_median(k) = ...
        PTT_median;

    Result_Rnumber(k) = ...
        length(R_loc);

    Result_PTTnumber(k) = ...
        length(PTT_clean);

    Result_NotchNumber(k) = ...
        NotchNumber;

    Result_NotchRate(k) = ...
        NotchRate;

end


%% =========================================================
% 18. 所有檔案比較
%% =========================================================

MeasurePosition = ...
    erase(Result_File, [".txt", "(1)"]);


Summary = table( ...
    MeasurePosition, ...
    Result_fs, ...
    Result_Duration, ...
    Result_HR, ...
    Result_PulseRate, ...
    Result_PTT_mean, ...
    Result_PTT_std, ...
    Result_PTT_median, ...
    Result_Rnumber, ...
    Result_PTTnumber, ...
    Result_NotchNumber, ...
    Result_NotchRate);


Summary.Properties.VariableNames = { ...
    'MeasurePosition', ...
    'Fs_Hz', ...
    'Duration_s', ...
    'ECG_HeartRate_bpm', ...
    'PPG_PulseRate_bpm', ...
    'PTT_Mean_ms', ...
    'PTT_SD_ms', ...
    'PTT_Median_ms', ...
    'Rwave_Number', ...
    'Valid_PTT_Number', ...
    'DicroticNotch_Number', ...
    'DicroticNotch_Rate_percent'};


disp(' ');
disp('============================================');
disp('ECG & PPG Analysis Result');
disp('============================================');

disp(Summary);
