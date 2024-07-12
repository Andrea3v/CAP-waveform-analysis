%% script to analyze and plot CAP. 

% Script works on .mat exported data from Patchmaster v>=2.9x
% data can be exported from both channels (V-mon on channel 2). 
% if only one channel, *A/D is found, value of pre-amplification used is required. 

% the script uses .1 - .9 ms as range for waveform offset & automatically detects the range for integration of the CAP area by: 
% 1) calculating the mean of the baseline waveforms
% 2a) using the prompt 'end of stim' [ms] to detect the beginning of the evoked CAP on the average waveforms
% 2b) if the rise of the evoked CAP response cannot be found, the user is required to add it manually by inserting a datatip to determine the left range for integration
% 3) using the peak latency of the baseline multiplied by a factor, e.g. 125%, to set the right end of the integration 

% range for calc of peak amplitude and latency has to be inserted manually by adding two datatips that determine the range where to find the waveform MAX

% the script will output the analyzed data as .xlsx in the subfolder 'Analysis' along with the final figure and the variables used for analysis
% The script won't run on exported data which might contain traces with both,one AND two channels, hence export only one type of data at time

clear all; clc 
[FileName,PathName] = uigetfile('*.mat','Select the MATLAB code file');
[~,NameNoExt,Ext] = fileparts(FileName);
ScriptFolder = cd(PathName);
S = load(fullfile(PathName,FileName));
off1 = 0.1;
off2 = 0.9;
Pix_SS = get(0,'screensize');
fig_ysize = 0.8; 
fig1=figure('Name',NameNoExt,...
    'Color', [0.99 0.99 0.99],...
    'Position', [0 round((Pix_SS(1,4)*(1-fig_ysize))) Pix_SS(1,3) round((Pix_SS(1,4)*fig_ysize))],...
    'Units', 'pixels',...
    'NumberTitle','off');
title('CAP');
%% separate channels
Ch1 = struct;
Ch2 = struct;
for f = fieldnames(S)'
    if str2num(f{1}(end)) == 1
        Ch1.(f{1}) = S.(f{1});
    elseif str2num(f{1}(end)) == 2
        Ch2.(f{1}) = S.(f{1});
    end
end
f = fieldnames(S)';
t0 = S.(char(f(1,1)))(1,1);
check = fieldnames(Ch2);
TF = double(isempty(check));
switch TF
    case 1
        prompt={'AD-channel! enter value for signal amplification used:',...
            'end of stimulation (ms)'...
            'Select baseline number of traces:',...
            'Enter right cut-off range for integration (n x peak latency):'};
        name = 'Amplification'; 
        defaultans = {'200','1.10','80:120','1.25'}; 
        answ = inputdlg(prompt,name,[1 40],defaultans);
        Mag = sprintf('%s*', answ{1,1});
        Mag = sscanf(Mag, '%f*');
        stims = str2num(answ{2,1});
        base = str2num(answ{3,1});
        rx = str2num(answ{4,1});
    case 0
        prompt={'end of stimulation (ms)'...
            'Select baseline number of traces:',...
            'Enter right cut-off range for integration (n x peak latency):'};
        name = 'Amplification'; 
        defaultans = {'1.05','1:120','1.25'}; 
        answ = inputdlg(prompt,name,[1 40],defaultans);
        stims = str2num(answ{1,1});
        base = str2num(answ{2,1});
        rx = str2num(answ{3,1});
        
end
%% channel 1
h = 1;
switch TF
    case 0
        s1 = subplot(1,2,1);
    case 1
end
for g = fieldnames(Ch1)'
    waveform = Ch1.(g{1});
    cap_time_ch1(h)=(waveform(1,1)-t0)/60;
    x = (waveform(:,1)-waveform(1,1))*1.00e+3;
    y = waveform(:,2)-mean(waveform(min(find(x> off1&x<off2)):max(find(x>off1&x<off2)),2));
    switch TF
        case 0
            y = y.*1.00e+6;
        case 1
            y = y.*(1000/Mag);        
    end
    Ch_1.(g{1}) = [x,y];
    y = sgolayfilt(y,4,23);
    plot(x,y);
    hold on;
    h = h+1;
end
xlim([0 10]);
switch TF
    case 0
        xlabel(s1,'ms');
        ylabel(s1,'nA');
        title('Current');
    case 1
        xlabel('ms');
        ylabel('/muV');
        title('Voltage');
end
%% channel 2
switch TF
    case 0
        h = 1;
        for g = fieldnames(Ch2)'
            s2 = subplot(1,2,2);
                waveform = Ch2.(g{1});
                cap_time_ch2(h)=(waveform(1,1)-t0)/60;
                x = (waveform(:,1)-waveform(1,1))*1.00e+3;
                y = waveform(:,2)-mean(waveform(min(find(x> off1&x<off2)):max(find(x>off1&x<off2)),2));
                Ch_2.(g{1}) = [x,y];
                y = sgolayfilt(y,4,23);
                plot(x,y);
                hold on
                h = h+1;
        end
        xlabel(s2,'ms'); 
        ylabel(s2,'mV');
        title('Voltage');
    case 1
end
xlim([0 10]);

%% grab range for CAP Area analysis
f_names = fieldnames(Ch_1)';
k = 1;
for j = base
    wav = Ch_1.(f_names{j});
    W(:,:,k) = wav;
    k = k + 1;
end
awb = mean(W,3); % average waveform baseline
x0 = awb(:,1);
y0 = awb(:,2);
try 
    logic0 = x0>stims;
awbCut = x0(logic0);
[p l] = findpeaks(-y0(logic0),'MinPeakProminence',.2);
%     plot(awb_x,awb_y)
    hold on
    plot(awbCut(l(1,1)),-p(1,1),'k^','MarkerFaceColor','g');
    left = awbCut(l(1,1));
catch
datacursormode on;
dcm_obj = datacursormode(gcf); 
set(dcm_obj,'Enable','on'); 
disp('select LEFT cut-off for CAP Area Analysis, then press ENTER')
ff = msgbox('select LEFT cut-off for CAP Area Analysis, then press ENTER','insert ONE datatip','help');
pause;
c_info2 = getCursorInfo(dcm_obj);
left = c_info2(1).Position(1,1);
end

    wav_x = x0(x0>left);
    wav_y = y0(x0>left);
    MMM = wav_x(wav_y == max(wav_y));
    latBase = mean(MMM);
right = rx*(mean(latBase));
hold on
rr = x0>=left&x0<right;
a= area(x0(rr),y0(rr),'FaceColor',[0 0 0],'LineStyle','none');
a.FaceAlpha = .3;

%% channel 1
h = 1;
for g = fieldnames(Ch_1)' 
    i = Ch_1.(g{1}); 
    x1=i(:,1);
    y1=i(:,2);
    cap_range=x1(x1>=left&x1<=right);
    cap_area_ch1(h)=trapz(cap_range,y1(x1>=left&x1<=right));
    h = h+1;
end
%% channel 2
switch TF
    case 0
        h = 1;
        for g = fieldnames(Ch_2)' 
            i = Ch_2.(g{1}); 
            x1=i(:,1);
            y1=i(:,2);
            cap_range=x1(x1>=left&x1<=right);
            cap_area_ch2(h)=trapz(cap_range,y1(x1>=left&x1<=right));
            h = h+1;
        end
    case 1
end
%% grab range for CAP amplitude and latency analysis
dcm_obj = datacursormode(gcf); 
set(dcm_obj,'Enable','on'); 
disp('select RANGE for CAP amplitude-latency Analysis, then press ENTER')
ff = msgbox('select RANGE for CAP amplitude-latency Analysis, then press ENTER','insert TWO datatips','help');
pause;
c_info3 = getCursorInfo(dcm_obj);
range0 = sort(vertcat(c_info3(1).Position(1,1),c_info3(2).Position(1,1)),'ascend')';
left0 = range0(1,1);
right0 = range0(1,2);
    if c_info3(1).Position(1,1)>c_info3(2).Position(1,1);
        xcap0 = c_info3(1).Target.XData;
        ycap0 = c_info3(1).Target.YData;
    else
        xcap0 = c_info3(2).Target.XData;
        ycap0 = c_info3(2).Target.YData;
    end
h = 1;
%% channel 1
for f = fieldnames(Ch_1)'
    i = Ch_1.(f{1}); 
    x2=i(:,1); 
    y2=i(:,2);
    cap_range=x2(x2>=left0&x2<=right0);
    cap_amplitude_ch1(h)=max(y2(x2>=left0&x2<=right0));
    cap_latency_ch1(h)=mean(x2(y2==max(y2(x2>=left0&x2<=right0))))-stims;
    h = h+1;
end
%% channel 2
switch TF
    case 0
        h = 1;
        for f = fieldnames(Ch_2)'
            i = Ch_2.(f{1}); 
            x2=i(:,1); 
            y2=i(:,2);
            cap_range=x2(x2>=left0&x2<=right0);
            cap_amplitude_ch2(h)=max(y2(x2>=left0&x2<=right0));
            cap_latency_ch2(h)=mean(x2(y2==max(y2(x2>=left0&x2<=right0))))-stims;
            h = h+1;
        end
    case 1
end

%% get all the var, plot and store analyzed data
cap_time_ch1 = cap_time_ch1';
cap_area_ch1 = cap_area_ch1';
cap_amplitude_ch1 = cap_amplitude_ch1';
cap_latency_ch1 = cap_latency_ch1';
switch TF 
    case 0
        cap_time_ch2 = cap_time_ch2';
        cap_area_ch2 = cap_area_ch2';
        cap_amplitude_ch2 = cap_amplitude_ch2';
        cap_latency_ch2 = cap_latency_ch2';
    case 1
end

fig2=figure('Name',NameNoExt,...
    'Color', [0.99 0.99 0.99],...
    'Position', [0 round((Pix_SS(1,4)*(1-fig_ysize)/2)) Pix_SS(1,3) round((Pix_SS(1,4)*fig_ysize))],...
    'Units', 'pixels',...
    'NumberTitle','off');
title('CAP Analysis');
hold on
dot_operator = char(8901);
%% subplot 1
switch TF 
    case 0
        sub1 = subplot(1,2,1);
    case 1
end
        sz =75;
        yyaxis left
        scatter(cap_time_ch1,cap_area_ch1,sz,'filled','p',...
            'LineWidth',.1,...
            'MarkerEdgeColor',[1 0 0],...
            'MarkerFaceColor',[.9 .3 .3]);
        hold on
switch TF 
    case 0
        sub1 = subplot(1,2,1);
    case 1
end
        sz =50;        
        yyaxis left
        scatter(cap_time_ch1,cap_amplitude_ch1,sz,'filled','o',...
            'LineWidth',.05,...
            'MarkerEdgeColor',[0, 0.4470, 0.7410],...
            'MarkerFaceColor',[0, 0.4470, 0.7410]);
        hold on
switch TF 
    case 0
        sub1 = subplot(1,2,1);
        xlabel(sub1,'min');
        ylabel(sub1,strcat('nA',dot_operator,'ms     ',' - ','     nA    '));
        title('Current')
    case 1
        xlabel('min');
        ylabel(strcat('mV',dot_operator,'ms     ',' - ','     /muV    ',' - ','   ms'));
        title('Voltage')
end
        sz =50;
        yyaxis right
        ylabel('ms (peak latency)')
        scatter(cap_time_ch1,cap_latency_ch1,sz,'filled','s',...
            'LineWidth',.05,...
            'MarkerEdgeColor',[0.850, 0.325, 0.098],...
            'MarkerFaceColor',[0.850, 0.325, 0.098]);
        hold on
%% subplot 2
switch TF 
    case 0
        sub2 = subplot(1,2,2);
        sz =100;
        yyaxis left;
        scatter(cap_time_ch2,cap_area_ch2,sz,'filled','p',...
            'LineWidth',.1,...
            'MarkerEdgeColor',[1 0 0],...
            'MarkerFaceColor',[.9 .3 .3]);
        sz =50;
        subplot(1,2,2);
        hold on
        yyaxis left
        xlabel(sub2,'min');
        ylabel(sub2,strcat('mV',dot_operator,'ms     ',' - ','     mV    '));
        scatter(cap_time_ch2,cap_amplitude_ch2,sz,'filled','o',...
            'LineWidth',.05,...
            'MarkerEdgeColor',[0, 0.4470, 0.7410],...
            'MarkerFaceColor',[0, 0.4470, 0.7410]);
        sz =50;
        subplot(1,2,2);
        yyaxis right;
        ylabel(sub2,'ms (peak latency)');
        scatter(cap_time_ch2,cap_latency_ch2,sz,'filled','s',...
            'LineWidth',.05,...
            'MarkerEdgeColor',[0.850, 0.325, 0.098],...
            'MarkerFaceColor',[0.850, 0.325, 0.098]);
        title('Voltage')
    case 1
end

legend('CAP Area','CAP Amplitude','CAP Latency','Location','best')

subdir = [PathName,'Analysis_files'];
    [status, msg, msgID] = mkdir(subdir);
    cd(subdir);
    warning('off',msgID);
file = strcat(subdir,'\',strtok(FileName,'.'),'.xlsx');
strange = strcat(num2str(range(1,1)),'_ms_',num2str(range(1,2)),'_ms');
strange0 = strcat(num2str(range0(1,1)),'_ms_',num2str(range0(1,2)),'_ms');
switch TF 
    case 0
        head = {'time_zero_I',strcat('cap_area_I_',strange),strcat('cap_amplitude_I_',strange0),strcat('cap_latency_I_',strange0),...
        'time_zero_V',strcat('cap_area_V_',strange),strcat('cap_amplitude_V_',strange0),strcat('cap_latency_V_',strange0)};
        T = table(cap_time_ch1,cap_area_ch1,cap_amplitude_ch1,cap_latency_ch1,...
          cap_time_ch2,cap_area_ch2,cap_amplitude_ch2,cap_latency_ch2);
        save([NameNoExt,'_Vars'],'Ch_1','Ch_2', 'cap_time_ch1', 'cap_area_ch1', 'cap_amplitude_ch1', 'cap_latency_ch1', ...
        'cap_area_ch1', 'cap_amplitude_ch1', 'cap_latency_ch1', 'strange', 'strange0');
    case 1
        head = {'time_zero_V',strcat('cap_area_V_',strange),strcat('cap_amplitude_V_',strange0),strcat('cap_latency_V_',strange0)};
        T = table(cap_time_ch1,cap_area_ch1,cap_amplitude_ch1,cap_latency_ch1);
        save([NameNoExt,'_Vars'],'Ch_1', 'cap_time_ch1', 'cap_area_ch1', 'cap_amplitude_ch1', 'cap_latency_ch1', 'strange', 'strange0');
end
writetable(T,file,'Sheet',1);
xlswrite(file,head,1,'A1:H1');
savefig(fig2,NameNoExt);
cd(ScriptFolder);
clearvars -except Ch_1 Ch_2 cap_time_ch1 cap_area_ch1 cap_amplitude_ch1 cap_latency_ch1 ...
     cap_time_ch2 cap_area_ch2 cap_amplitude_ch2 cap_latency_ch2 NameNoExt PathName strange strange0