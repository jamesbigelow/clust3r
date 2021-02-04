function [ cidx_out, waves_out, spiketimes_out, sortinfo] = clust3r( waves0, spiketimes0, opts, cidx0 )
% NAME: CLUST3R
%
% DESCRIPTION: Interactive manual spike sorting program for extracellular 
% spike waveform data. Clusters are defined by freehand contours drawn by 
% user in 3D feature space. Semi-automated sorting is possible via 
% Template matching functions.
%
% INPUT:
% waves: (required) Spike waveforms, N x K array (spike index x n samples)
% spiketimes: (required) Spiketimes, N length vector (default sec; see opts)
% opts: (optional) options structure with details about user and data
% cidx0: (optional) cluster indices for loading sorted data
%
% OUTPUT:
% cidx: Cluster index, each cell contains indices of clustered spikes. Up
% to 5 clusters may be defined. If no clusters are defined, or if user
% exits without clicking 'Exit & Save' button, cidx is an empty cell array.
% waves_out: waveforms retained during sorting
% spiketimes_out: spiketimes retained during sorting
% sortinfo: structure with details about sorting
%
% INTERNAL FUNCTIONS: 
% plott3r, wav3r, hist3r, fsett3r (slid3r, slic3r), hotk3ys, alldelet3r
% add3r (dragg3r), remov3r (dragg3r), combin3r (check3r), allmatch3r,
% pidx3r (check3r), allpidx3r, xcorr3r, resett3r, sav3r, solopidx3r,
% delet3r, match3r, origin3r, explor3r (click3r), pcview3r, view3r, norm3r,
% spkisi3r, featur3r, figur3r, stat3r, color3r
%
% EXTERNAL FUNCTION: 
% clust3r_stats
%  
% (C) James Bigelow (jamesbigelow at gmail dot com)

% Generate options structure if not provided ---------
if ~exist('opts','var')
    opts.sorterID = 'AN';
    opts.experimentTime = '2015-10-21_block-1';
    opts.savepath = pwd; % Where to save data files 
    opts.figpath = opts.savepath; % Where to save figures
    opts.channel = 0;
    opts.nWavSamples = size(waves0,2); % number of samples per channel (for multichannel waveform data)
    opts.tspkUnits = 's'; % 's' for seconds, 'ms' for milliseconds, 'samples' for samples
    opts.fs = ''; % sample rate (for documentation purpose only)
end
opts.version = '2018.08.23';
sort_time = datestr(now,'yyyy-mm-dd_HH-MM-ss');

% Initialize persistent internal variables ---------
if ~iscolumn(spiketimes0); spiketimes0 = spiketimes0'; end
persistent waves spiketimes uu vv PC1_ PC2_ PC3_ PC1 PC2 PC3 Peak Valley PeakValley Variability Spiketimes PC1_0 PC2_0 PC3_0 Peak_0 Valley_0 PeakValley_0 Variability_0 Spiketimes_0 nspk cidx pidx Xid Yid Zid X Y Z X_0 Y_0 Z_0 fstrings rstrings
persistent Slice1 Slice2 Slice3 Slice1_0 Slice2_0 Slice3_0 Template1 Template2 Template3 Template1_0 Template2_0 Template3_0
persistent fign chan figpath cluster_statistics; fign = 0; chan = opts.channel; figpath = opts.figpath; cluster_statistics = [];
waves = waves0; spiketimes = spiketimes0; % Store originals
fstrings = {'PC1','PC2','PC3','Spiketimes','PeakValley','Valley','Peak','Variability','Slice1','Slice2','Slice3','Template1','Template2','Template3'}; %  Added fstrings to facilitate updating feature popup strings
rstrings = []; %  Added runstrings variable to enable post sorting command evaluation 
persistent existing_sort_flag cidx_existing %  added variables to facilitate viewing vs. editing existing sorts (see sav3r) 

% GUI ==================================================
close all
saveflag = false;
nsmpl = 200; % Sample of waveforms to plot
colorvec = [.4,.4,.4; 1,0,0; 0,0,1; 1,0,1; 0,1,1; 0,1,0;]; %  Changed cluster 1 color to lighter gray for better distinction from cluster 3
plotcolor = 'none'; %  added for controlling background color of plots
fntsz = 12; if max(get(groot,'Screensize'))<2000; fntsz = fntsz/2; end
gf = figure('KeyPressFcn', @hotk3ys,'toolbar','none','menu','none','Tag','GUI','Position', get(0, 'Screensize')); % GUI figure;  added hotk3ys function, removed toolbar and menu to avoid problems with hotk3ys (see fix3r); Added GUI tag for preservation in sav3r function.
% Maximize 
try
    drawnow;
    set(gf, 'WindowStyle', 'normal');
    oldState = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jFig = get(handle(gf), 'JavaFrame');
    jFig.setMaximized(true);
catch
end
set(gf,'name',['clust3r (v ' opts.version ')'],'numbertitle','off','UserData',0);

T1=110;T2=45; % Default view % T1=-60;T2=60;
persistent ms tf msp ws hs origins
persistent wxh wyl wyh

% Control buttons ---------------------------------------------------------
sbx = .001; sbw = .06; bfsz=8; % Button font size  added sbx, sbw to facilitate changes to side button positions and sizes 
g1top = .92; bspace = .03; tspace = .015; %  added g1top, g2top, g3top, g4top, g5top, bspace, tspace to facilitate changes to side button positions
% Group 1
sb(1) = uicontrol('Style','text','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g1top+.01 sbw .02],'String','Feature 1');
sb(2) = uicontrol('Style','popup','UserData',1,'String',fstrings,... %  Updated with fstrings variable for facility in updating feature popup strings [sb(2), sb(4), sb(6)]
    'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position', [sbx g1top sbw .02],'Callback',@fsett3r);
sb(3) = uicontrol('Style','text','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g1top+.01-bspace sbw .02],'String','Feature 2');
sb(4) = uicontrol('Style','popup','UserData',2,'String',fstrings,...
    'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position', [sbx g1top-bspace sbw .02],'Callback',@fsett3r);
sb(5) = uicontrol('Style','text','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g1top+.01-bspace*2 sbw .02],'String','Feature 3');
sb(6) = uicontrol('Style','popup','UserData',3,'String',fstrings,...
    'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position', [sbx g1top-bspace*2 sbw .02],'Callback',@fsett3r);
sb(7) = uicontrol('Style','togglebutton','String','3D origins',...
    'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g1top-bspace*3 sbw .02],'Callback',@origin3r);
sb(27) = uicontrol('Style','togglebutton','String','White background',...
    'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g1top-bspace*4 sbw .02],'Callback',@color3r); %  added button to control plot background color
% Group 2
g2top = .75;
sb(8) = uicontrol('String','[A]dd cluster','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g2top sbw .02],'Callback',@add3r); %  side button mods to string and size/position to indicate hotkeys
sb(9) = uicontrol('String','[R]emove waves','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g2top-bspace sbw .02],'Callback',@remov3r);
sb(10) = uicontrol('String','Explore waves','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g2top-bspace*2 sbw .02],'Callback',@explor3r);
sb(11) = uicontrol('String','[E]xcise waves','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g2top-bspace*3 sbw .02],'Callback',@excis3r);
% Group 3
g3top = .60;
sb(12) = uicontrol('String','Combine clusters','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top sbw .02],'Callback',@combin3r);
sb(13) = uicontrol('String','Plot multiple','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace sbw .02],'Callback',@pidx3r);
sb(14) = uicontrol('String','[P]lot all','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*2 sbw .02],'Callback',@allpidx3r);
sb(15) = uicontrol('String','[D]elete all','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*3 sbw .02],'Callback',@alldelet3r);
sb(16) = uicontrol('String','[T]emplate match all','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*4 sbw .02],'Callback',@allmatch3r);
sb(17) = uicontrol('String','[X]Cross correlations','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*5 sbw .02],'Callback',@xcorr3r);
sb(18) = uicontrol('String','[C]luster statistics','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*6 sbw .02],'Callback',@stat3r);
sb(19) = uicontrol('String','View PC','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*7 sbw .02],'Callback',@pcview3r);
sb(20) = uicontrol('String','Redefine channels','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g3top-bspace*8 sbw .02],'Callback',@channel3r);
% Group 4
g4top = .30;
sb(21) = uicontrol('String','Save [f]igure','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g4top sbw .02],'Callback',@figur3r);
sb(22) = uicontrol('String','Reset waves/features','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g4top-bspace sbw .02],'Callback',@resett3r);
sb(23) = uicontrol('String','[Enter] Save & Exit','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g4top-bspace*2 sbw .02],'Callback',@sav3r);
% Group 5
g5top = .19;
sb(24) = uicontrol('Style','text','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g5top+.015 .05 .02],'String','Sorter ID');
sb(25) = uicontrol('Style','edit','String',opts.sorterID,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position', [sbx g5top sbw .02]);

sb(26) = uicontrol('String','[Q]uit w/o saving','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g1top+bspace*1.5 sbw .02],'Callback',@exit3r); %  added button to quit w/o saving, deliberately distant from Save/proceed button 

% Group 6  
%{
g6top = .16;
sb(28) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top .075 .02],'String','UNIT SORT CODES');
sb(38) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace .075 .02],'String','[Enter in colored boxes]');
sb(29) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*2 .075 .02],'String','00 - Raw multiunit');
sb(30) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*3 .075 .02],'String','01 - Clean multiunit');
sb(31) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*4 .075 .02],'String','11 - Excellent cluster #1');
sb(32) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*5 .075 .02],'String','12 - Excellent cluster #2');
sb(33) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*6 .075 .02],'String','21 - Decent cluster #1');
sb(34) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*7 .075 .02],'String','22 - Decent cluster #2');
sb(35) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*8 .075 .02],'String','31 - Poor cluster #1');
sb(36) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*9 .075 .02],'String','32 - Poor cluster #2');
sb(37) = uicontrol('Style','text','HorizontalAlignment', 'left','Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[sbx g6top-tspace*10 .075 .02],'String','99 - Artifact');
%}

% 2D view buttons ---------------------------------------------------------
tb(1) = uicontrol('String','[F1] Enable view','UserData',1,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[.93 .955 .05 .02],'Callback',@view3r);
tb(2) = uicontrol('String','[F2] Enable view','UserData',2,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[.93 .706 .05 .02],'Callback',@view3r);
tb(3) = uicontrol('String','[F3] Enable view','UserData',3,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[.93 .457 .05 .02],'Callback',@view3r);
tb(4) = uicontrol('String','[F4] Enable view','UserData',4,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[.93 .209 .05 .02],'Callback',@view3r);

% Cluster buttons ------------------------------------------------------------
uy=.2; uw=.025; uh=.025;
for uii = 1:6
    ux = .08+(uii*.12)-.12;
    ub1(uii) = uicontrol('String',['Select [' num2str(uii) ']'],'UserData',uii,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[ux uy uw*1.4 uh],'Callback',@solopidx3r); %  mod to string and size/position to indicate hotkeys
    ub2(uii) = uicontrol('String','Delete','UserData',uii,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[ux+.036 uy uw*.9 uh],'Callback',@delet3r); %.035
    ub3(uii) = uicontrol('String','Match','UserData',uii,'Units','normalized','Fontsize',bfsz,'FontUnits','normalized','Position',[ux+.06 uy uw*.9 uh],'Callback',@match3r); %.07
    ub4(uii) = uicontrol('Style','edit','String','','UserData',uii,'Units','normalized','Fontsize',ceil(bfsz * 1.2),'FontUnits','normalized','Position',[ux+.0825 uy uw*.72 uh * 1.0], 'BackgroundColor', colorvec(uii, :) ./ 2 + [.3 .3 .3]);
end
if ~isfield(opts,'qual')
    % set(ub4(1),'String','00'); %  Updated default to 00 for 'uninspected multiunit'
else
    for uii = 1: length( opts.qual )
        set(ub4(uii),'String',opts.qual{uii}); %  Updated to populate unit designation codes from previously-sorted data
    end
end

% Plot unsorted data ==================================================
featur3r;plott3r;plott3r2D;wav3r;hist3r;

% Wait until user clicks Exit or Save & Exit, update cidx and close -------
if isfield(opts,'autosave') %  Added to enable autosave functionality (see external function time_limit_sorted_folder)
    sav3r; saveflag = true; 
else
    waitfor(gf,'UserData',1);
end
plott3r;plott3r2D;wav3r;hist3r;

if saveflag

    % Overwrite previous sorted data for this channel, if any exists
    %  First, delete data files in sorted folder...
    fids = dir( fullfile( opts.savepath, sprintf('*ch%02d*.mat',opts.channel ))); % updated to reflect file naming behavior
    for ff = 1:length(fids)
        if fids(ff).date < datetime( sort_time, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss' )
            delete(fullfile(opts.savepath, fids(ff).name));
        end
    end
    %  Second, delete old figures in plots folder...
    fids = dir( fullfile( opts.figpath, sprintf('*Channel*%02d*.png',opts.channel )));
    for ff = 1:length(fids)
        if fids(ff).date < datetime( sort_time, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss' )
            delete(fullfile(opts.figpath, fids(ff).name));
        end
    end
    figur3r([],[], true); % replaced figur3r call here so that if saveflag is true, the plot for the new/edited sort is generated after the old plot is deleted by the lines above
    
    opts.sorterID = get ( sb(25) ,'string') ; %  store user ID if updated during sorting
    
    % Store cluster statistics and sort options in 'sortinfo' structure
    cidx_out = cidx(2:max(find(~cellfun(@isempty,cidx)))); % Index only sorted clusters
    waves_out = waves;
    spiketimes_out = spiketimes;
    % [cluster_stats_f] = cluster_stats(cidx(2:max(find(~cellfun(@isempty,cidx)))),[Y Z],waves,spiketimes,0,[],colorvec);
    [cluster_stats_f] = clust3r_stats(cidx(2:max(find(~cellfun(@isempty,cidx)))),[Y Z],waves,spiketimes,0,[],colorvec);
    for ii = 1:length(cidx)
        if ~isempty(cidx{ii})
            cluster_stats_f.qual{ii} = get(ub4(ii),'String');            
        end
    end
    cluster_stats_f.features{1}= Yid;
    cluster_stats_f.features{2}= Zid;
    if isempty(cluster_statistics)
        cluster_statistics = cluster_stats_f;
    else
        cluster_statistics(end+1) = cluster_stats_f;
    end
    sortinfo.cluster_statistics = cluster_statistics;
    sortinfo.opts = opts;
    sortinfo.sort_time = sort_time;
    
    % Save mat file for each cluster containing 'nev' structure with
    % spiketimes, waves, sort indices and 'sortinfo' structure
    for u = 1:length(cidx)
        if ~isempty(cidx{u})
            
            nev.spiketimes = spiketimes(cidx{u});
            nev.waves = waves(cidx{u},:);
            nev.sortindices = cidx{u};
            nev.sortinfo = sortinfo;
            filename = ['sorted_' opts.experimentTime '_ch' sprintf('%02d',opts.channel) '_cluster' sprintf('%02d',u) '.mat']; %  Updated file naming behavior per group discussion
            
            save( fullfile(opts.savepath, filename),'nev');
        end
    end
    disp(['Channel ' num2str(opts.channel) ' Done'])
else
    
    cidx_out = {};
    waves_out = [];
    spiketimes_out = [];
    if isempty(rstrings) % Added rstrings to sortinfo to enable post sorting command evaluation
        sortinfo = [];
    else
        sortinfo.rstrings = rstrings;
    end
    
end
close all

% Initialize wave features ==================================================

% FEATUR3R
    function featur3r
        % Waveform PCA and other features -----------------------------------------
        [uu,~,vv] = svd(waves,'econ');
        % Force to first quadrant
        if mean(uu(:,1))<0; uu(:,1)=-1*uu(:,1); vv(:,1)=-1*vv(:,1); end
        if mean(uu(:,2))<0; uu(:,2)=-1*uu(:,2); vv(:,2)=-1*vv(:,2); end
        if mean(uu(:,3))<0; uu(:,3)=-1*uu(:,3); vv(:,3)=-1*vv(:,3); end
        PC1_ = uu(:,1); PC2_ = uu(:,2); PC3_ = uu(:,3);
        Peak = max(waves,[],2); Valley = min(waves,[],2); PeakValley = Peak-Valley;
        % Variability = mean(diff(waves,1,2),2);
        for w = 1:length(waves)
            D = diff(waves(w,:))>0;
            D0 = xor(D(1:end-1),D(2:end));
            D00 = abs(diff(waves(w,find(D0))));
            Variability(w,1) = 1/( max(D00) / median(D00)) ;
        end
        % Normalize features between 0 and 1 for plot3 facility
        PC1 = norm3r(PC1_); PC2 = norm3r(PC2_); PC3 = norm3r(PC3_);
        Spiketimes = norm3r(spiketimes);
        Peak = norm3r(Peak); Valley = norm3r(Valley); PeakValley = norm3r(PeakValley); Variability = norm3r(Variability);
        % Normalized 0 values for optional origin lines
        PC1_0 = norm3r([0; uu(:,1)]); PC1_0 = PC1_0(1);
        PC2_0 = norm3r([0; uu(:,2)]); PC2_0 = PC2_0(1);
        PC3_0 = norm3r([0; uu(:,3)]); PC3_0 = PC3_0(1);
        Spiketimes_0 = min(Spiketimes);
        Peak_0 = min(Peak); Valley_0 = min(Valley); PeakValley_0 = min(PeakValley); Variability_0 = min(Variability);
        nspk = length(Spiketimes);
        
        if ~exist('cidx0','var')
            existing_sort_flag = false; %  added to facilitate reviewing process 
            cidx = cell(1,6); % Cluster index
            cidx{1} = 1:nspk; % Initialize cluster 1 with unsorted spikes
            pidx = 1; % Plot index
        else
            existing_sort_flag = true; %  added to facilitate reviewing process 
            cidx_existing = cidx0; %  added to check if anything was changed during review process
            if length(cidx0) == 6
                cidx = cidx0; clear cidx0 %  Clear after importing to avoid resett3r bug
                pidx = find(~cellfun('isempty',cidx));
            else
                cidx1 = cidx0; clear cidx0 %  Clear after importing to avoid resett3r bug
                cidx = cell(1,6); % Cluster index
                cidx{1} = setdiff(1:nspk,[cidx1{:}]);
                for cl = 1:length(cidx1)
                    cidx{cl+1} = cidx1{cl};
                end
                pidx = find(~cellfun('isempty',cidx));
            end
        end
        % Default 3D feature space
        Xid = 'Spiketimes';     X = eval(Xid);      X_0 = eval(strcat(Xid,'_0'));
        Yid = 'PC1';            Y = eval(Yid);      Y_0 = eval(strcat(Yid,'_0'));
        Zid = 'PC2';            Z = eval(Zid);      Z_0 = eval(strcat(Zid,'_0'));
        origins=1;
        wxh = size(waves,2); wyl = quantile(min(waves,[],2),.002); wyh = quantile(max(waves,[],2),.998);
        % Update popup menu strings (JB updated 2018-08-13)
        fidx2 = find( strcmp( fstrings, Xid) ); set( sb(2), 'Value', fidx2);
        fidx4 = find( strcmp( fstrings, Yid) ); set( sb(4), 'Value', fidx4);
        fidx6 = find( strcmp( fstrings, Zid) ); set( sb(6), 'Value', fidx6);
    end

    function plott3r2D
        % 2D feature space plots --------------------------------------------------
        try delete(tf); catch; end
        tw = .16; th = .2;
        msz2 = 2; if max(get(groot,'Screensize'))<2000; msz2 = msz2/10; end
        
        tf(1) = subplot('Position',[.82 .78 tw th]);
        hold on
        for i = 1:6
            if ~isempty(cidx{i})
                plot(uu(cidx{i},1),uu(cidx{i},2),'k.','markersize',msz2,'MarkerEdgeColor',colorvec(i,:));
            end
        end
        hold off
        line([min(xlim) max(xlim)],[0 0],'color','k','linestyle',':');
        line([0 0],[min(ylim) max(ylim)],'color','k','linestyle',':');
        xlabel('PC1','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized'); ylabel('PC2','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized');
        set(tf(1),'box','on','xtick',[],'ytick',[]); set(tf(1),'color',plotcolor); %  Updated to control plot background color
        
        tf(2) = subplot('Position',[.82 .53 tw th]);
        hold on
        for i = 1:6
            if ~isempty(cidx{i})
                plot(uu(cidx{i},1),PeakValley(cidx{i}),'k.','markersize',msz2,'MarkerEdgeColor',colorvec(i,:));
            end
        end
        hold off
        line([0 0],[min(ylim) max(ylim)],'color','k','linestyle',':');
        xlabel('PC1','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized'); ylabel('PeakValley','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized');
        set(tf(2),'box','on','xtick',[],'ytick',[]); set(tf(2),'color',plotcolor); %  Updated to control plot background color
        
        tf(3) = subplot('Position',[.82 .28 tw th]);
        hold on
        for i = 1:6
            if ~isempty(cidx{i})
                plot(Spiketimes(cidx{i}),uu(cidx{i},1),'k.','markersize',msz2,'MarkerEdgeColor',colorvec(i,:));
            end
        end
        hold off
        line([min(xlim) max(xlim)],[0 0],'color','k','linestyle',':');
        xlabel('Spiketimes','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized'); ylabel('PC1','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized');
        set(tf(3),'box','on','xtick',[],'ytick',[]); set(tf(3),'color',plotcolor); %  Updated to control plot background color
        
        tf(4) = subplot('Position',[.82 .03 tw th]);
        hold on
        for i = 1:6
            if ~isempty(cidx{i})
                plot(Spiketimes(cidx{i}),PeakValley(cidx{i}),'k.','markersize',msz2,'MarkerEdgeColor',colorvec(i,:));
            end
        end
        hold off
        xlabel('Spiketimes','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized'); ylabel('PeakValley','FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized');
        set(tf(4),'box','on','xtick',[],'ytick',[]); set(tf(4),'color',plotcolor); %  Updated to control plot background color
    end

% Dynamic plot functions ==================================================

% PLOTT3R
    function plott3r
        % 3D feature space plot -------------------------------------------
        try delete(ms); catch; end
        ms = subplot('Position',[.08 .26 .7 .72]);
        hold on
        for pii = 1:length(pidx)
            if pidx(pii)==1;msz=8;else msz=10;end
            if max(get(groot,'Screensize'))<3000; msz = msz/1.5; end
            if ~isempty(cidx{pidx(pii)})
                msp = plot3(X(cidx{pidx(pii)}),Y(cidx{pidx(pii)}),Z(cidx{pidx(pii)}),'.','markersize',msz,'MarkerEdgeColor',colorvec(pidx(pii),:));
            end
        end
        if origins
            lwo = 2; if max(get(groot,'Screensize'))<2000; lwo = 1; end
            line([0 1],[Y_0 Y_0],[Z_0 Z_0],'color','k','linewidth',lwo,'linestyle',':','HitTest','off');
            line([X_0 X_0],[0 1],[Z_0 Z_0],'color','k','linewidth',lwo,'linestyle',':','HitTest','off');
            line([X_0 X_0],[Y_0 Y_0],[0 1],'color','k','linewidth',lwo,'linestyle',':','HitTest','off');
        end
        hold off; view(T1,T2);
        xlabel(Xid, 'FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized');ylabel(Yid, 'FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized');zlabel(Zid, 'FontWeight','bold','Fontsize',fntsz,'FontUnits','normalized')
        set(ms,'box','on','xtick',[],'ytick',[],'ztick',[],'xlim',[-.05 1.05],'ylim',[-.05 1.05],'zlim',[-.05 1.05]); set(ms,'color',plotcolor); %  Updated to control plot background color
        title(['Channel ' num2str(opts.channel)],'fontsize',fntsz*.75)
        rotate3d(ms,'on');  %  Modified rotate3d behavior so that it applies only to the 3D figure; removed all other instances of 'rotate3d on'
        % Restore KeyPressFcn for hotk3ys following rotate3d on
        % https://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
        hManager = uigetmodemanager(gf);
        try
            set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
        catch
            [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
        end
        set(gf, 'WindowKeyPressFcn', []);
        set(gf, 'KeyPressFcn', @hotk3ys);
    end

% WAV3R
    function wav3r
        % Waveform plots --------------------------------------------------
        try delete(ws); catch; end
        wy = .11; ww = .1; wh = .09;
        for wii = 1:6
            wx = .08+(wii*.12)-.12;
            ws(wii) = subplot('Position',[wx wy ww wh]);
            if ~isempty(cidx{wii})
                hold on
                for wjj = 1:nsmpl
                    plot(waves(cidx{wii}(randi(length(cidx{wii}))),:),'color',[.6,.6,.6]);
                end
                plot(mean(waves(cidx{wii},:)),'color',colorvec(wii,:),'linewidth',3);
                if length(waves(cidx{wii})) > 30
                    plot(mean(waves(cidx{wii},:))+std(waves(cidx{wii},:)),'color',colorvec(wii,:),'linewidth',1.5,'linestyle',':');
                    plot(mean(waves(cidx{wii},:))-std(waves(cidx{wii},:)),'color',colorvec(wii,:),'linewidth',1.5,'linestyle',':');
                end
                wtext = strcat({'Cluster '},{num2str(wii)}); %  updated naming scheme designates everything as cluster
                text(max(xlim)*.75,wyh*.8,wtext,'Color',colorvec(wii,:),'fontweight','bold',...
                    'FontSize',fntsz,'FontUnits','normalized','interpreter','none','horizontalAlignment','right');
                hold off
            else
                plot([]);
            end
            set(ws(wii),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        end
    end

% HIST3R
    function hist3r
        try delete(hs); catch; end
        hy = .01; hw = .1; hh = .09;
        for hii = 1:6
            hx = .08+(hii*.12)-.12;
            hs(hii) = subplot('Position',[hx hy hw hh]);
            if numel(cidx{hii})>2 %  Updated to prevent crash with only one value
                [bins,isi1] = spkisi3r(spiketimes(sort(cidx{hii})));
                bar(bins,1,'FaceColor',colorvec(hii,:),'EdgeColor','none'); xlim([0.5 length(bins)+.5])
                line([18 18],[0 max(ylim)],'color','k','linestyle','--');
                text(max(xlim)*.95,max(ylim)*.9,strcat(['Short ISI: ',num2str(isi1,'%.1f')],'%'),...
                    'Color','k','fontweight','bold','FontSize',fntsz,'FontUnits','normalized','interpreter','none','horizontalAlignment','right');
                text(max(xlim)*.95,max(ylim)*.8,strcat(['n = ',num2str(length(cidx{hii}))]),...
                    'Color','k','fontweight','bold','FontSize',fntsz,'FontUnits','normalized','interpreter','none','horizontalAlignment','right');
                
                switch opts.tspkUnits
                    case 's'
                        spkHz = num2str( length(cidx{hii}) / ( (spiketimes(end)-spiketimes(1))/1000)/1000 ,'%.1f');
                    case 'ms'
                        spkHz = num2str( length(cidx{hii}) / ( (spiketimes(end)-spiketimes(1))/1000) ,'%.1f');
                    case 'samples'
                        spkHz = num2str( length(cidx{hii}) / ( (spiketimes(end)-spiketimes(1))/opts.fs)/1000 ,'%.1f');
                end
                text(max(xlim)*.95,max(ylim)*.7,['Hz: ' spkHz ],...
                    'Color','k','fontweight','bold','FontSize',fntsz,'FontUnits','normalized','interpreter','none','horizontalAlignment','right');
            else
                plot([]);
            end
            set(hs(hii),'color','none','box','on','xtick',[12 24 36 48],'ytick',[],'xticklabel',[],'TickLength',[0.025, 0.01],'Layer','top');
        end
    end

% Button callback functions ===============================================

% FSETT3R
    function fsett3r(source,~)
        fid = source.UserData;
        fval = source.Value;
        fstrings = source.String;
        persistent sf sld slca slcb rfig euc
        
        switch fid
            case 1
                set(sb(2), 'Enable', 'off'); drawnow; set(sb(2), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
                if strfind(fstrings{fval},'Slice')
                    slic3r
                    Xid = fstrings{fval}; X = eval(Xid); X_0 = eval(strcat(Xid,'_0'));
                elseif strfind(fstrings{fval},'Template')
                    templat3r
                    Xid = fstrings{fval}; X = eval(Xid); X_0 = eval(strcat(Xid,'_0'));
                else
                    Xid = fstrings{fval}; X = eval(Xid); X_0 = eval(strcat(Xid,'_0'));
                end
            case 2
                set(sb(4), 'Enable', 'off'); drawnow; set(sb(4), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
                if strfind(fstrings{fval},'Slice')
                    slic3r
                    Yid = fstrings{fval}; Y = eval(Yid); Y_0 = eval(strcat(Yid,'_0'));
                elseif strfind(fstrings{fval},'Template')
                    templat3r
                    Yid = fstrings{fval}; Y = eval(Yid); Y_0 = eval(strcat(Yid,'_0'));
                else
                    Yid = fstrings{fval}; Y = eval(Yid); Y_0 = eval(strcat(Yid,'_0'));
                end
            case 3
                set(sb(6), 'Enable', 'off'); drawnow; set(sb(6), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
                if strfind(fstrings{fval},'Slice')
                    slic3r
                    Zid = fstrings{fval}; Z = eval(Zid); Z_0 = eval(strcat(Zid,'_0'));
                elseif strfind(fstrings{fval},'Template')
                    templat3r
                    Zid = fstrings{fval}; Z = eval(Zid); Z_0 = eval(strcat(Zid,'_0'));
                else
                    Zid = fstrings{fval}; Z = eval(Zid); Z_0 = eval(strcat(Zid,'_0'));
                end
        end
        
        [T1,T2]=view(ms);plott3r; % Update plot
        
        function slic3r(varargin)
            % Get slice from slider GUI
            rotate3d off; pan off; zoom off; datacursormode off; brush off; plotedit off
            if numel(pidx)>1; spidx = max(pidx); else spidx = pidx; end
            sf = figure('units','normalized','position',[.25,.4,.7,.55],'toolbar','none','menu','none','name','','numbertitle','off');
            title('[CLOSE WINDOW TO IMPLEMENT SLICE CHOICE]'); % added instruction title to clarify user must close slice subplot window to implement choice and continue sorting 
            hold on
            for sii = 1:nsmpl
                plot(waves(cidx{spidx}(randi(length(cidx{spidx}))),:),'color',[.6,.6,.6]);
            end
            plot(mean(waves(cidx{spidx},:)),'color',colorvec(spidx,:),'linewidth',3);
            if length(waves(cidx{spidx})) > 30
                plot(mean(waves(cidx{spidx},:))+std(waves(cidx{spidx},:)),'color',colorvec(spidx,:),'linewidth',1.5,'linestyle',':');
                plot(mean(waves(cidx{spidx},:))-std(waves(cidx{spidx},:)),'color',colorvec(spidx,:),'linewidth',1.5,'linestyle',':');
            end
            set(sf.CurrentAxes,'color','none','box','off','xcolor','none','ycolor','none','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
            sld = uicontrol('Style','slider','Min',1,'Max',wxh,'Value',round(wxh/2),'units','normalized',...
                'Position', [.13 0 .775 .1],'Callback', @slid3r);
            slca = line([round(wxh/2) round(wxh/2)],[wyl wyh],'color','k','linewidth',2);
            set(sf,'CloseRequestFcn',@slid3r1);
            set(sf,'UserData',0)
            waitfor(sf,'UserData',1);
            % Update Slice data
            stmp = waves(:,slcb);
            stmp_0 = norm3r([0; stmp]); stmp_0 = stmp_0(1);
            assignin('caller', fstrings{fval}, norm3r(stmp));
            assignin('caller', strcat(fstrings{fval},'_0'), stmp_0);
        end
        function slid3r(source,~)
            delete(slca)
            slca = line([round(source.Value) round(source.Value)],[wyl wyh],'color','k','linewidth',2);
        end
        function slid3r1(varargin)
            set(sf,'CloseRequestFcn','');
            slcb = round(get(sld,'Value'));
            set(sf,'UserData',1);
            delete(sf);
        end
        function templat3r(varargin)
            if length(find(~cellfun(@isempty,cidx)))<2;errordlg('Add cluster to define template');end
            rotate3d off; pan off; zoom off; datacursormode off; brush off; plotedit off
            rfw=200; rfh=200;
            rfig = figure('units','pixels','position',[100,100,rfw,rfh],'toolbar','none','menu','none');
            set(rfig,'name','Which template?','numbertitle','off');
            
            % Template radio buttons
            rrb(1) = uicontrol('style','radiobutton','units','pixels',...
                'position',[10,rfh*.86,70,15],'string','Unit 1');
            rrb(2) = uicontrol('style','radiobutton','units','pixels',...
                'position',[10,rfh*.70,70,15],'string','Unit 2');
            rrb(3) = uicontrol('style','radiobutton','units','pixels',...
                'position',[10,rfh*.54,70,15],'string','Unit 3');
            rrb(4) = uicontrol('style','radiobutton','units','pixels',...
                'position',[10,rfh*.38,70,15],'string','Unit 4');
            rrb(5) = uicontrol('style','radiobutton','units','pixels',...
                'position',[10,rfh*.22,70,15],'string','Unit 5');
            
            % Select button
            uicontrol('style','pushbutton','units','pixels',...
                'position',[100,rfh*.45,70,20],'string','Select',...
                'callback',@select3r);
            set(rfig,'UserData',0)
            waitfor(rfig,'UserData',1);
            euc_0 = norm3r([0; euc]); euc_0 = euc_0(1);
            assignin('caller', fstrings{fval}, norm3r(euc));
            assignin('caller', strcat(fstrings{fval},'_0'), euc_0);
            function select3r(varargin)
                rrbs = get(rrb,'Value');
                selected = find([rrbs{:}]);
                if ~isempty(selected)
                    tmplt = mean(waves(cidx{selected+1},:));
                    euc = nan(length(waves),1);
                    for tjj = 1:length(waves)
                        euc(tjj) = 1/log( ( sqrt(sum((waves(tjj,:)-tmplt).^2)) ) );
                    end
                    set(rfig,'UserData',1);
                    delete(rfig);
                else
                    set(rfig,'UserData',1);
                    delete(rfig);
                end
            end
        end
    end

% PCVIEW3R
    function pcview3r(source,~)
        set(sb(19), 'Enable', 'off'); drawnow; set(sb(19), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        % Quantiles
        PCahi = quantile(uu(:,1),.95); PCalo = quantile(uu(:,1),.05);
        PCbhi = quantile(uu(:,2),.95); PCblo = quantile(uu(:,2),.05);
        PCchi = quantile(uu(:,3),.95); PCclo = quantile(uu(:,3),.05);
        % Wave samples
        wsPCahi = waves(find(uu(:,1)>PCahi),:); wsPCalo = waves(find(uu(:,1)<PCalo),:);
        wsPCbhi = waves(find(uu(:,2)>PCbhi),:); wsPCblo = waves(find(uu(:,2)<PCblo),:);
        wsPCchi = waves(find(uu(:,3)>PCchi),:); wsPCclo = waves(find(uu(:,3)<PCclo),:);
        % Y limits
        pyl = min(min(vv(:,1:3)))*1.1; pyh = max(max(vv(:,1:3)))*1.1;
        % Plot
        pfig = figure('Units','normalized','position',[.2,.2,.45,.6],'toolbar','none','menu','none','name','','numbertitle','off'); %  Added handle to close 
        % PC subplots
        ps(1) = subplot(331);
        plot(vv(:,1),'k','linewidth',1.5);title('PC1')
        set(ps(1),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[pyl pyh]);
        ps(2) = subplot(332);
        plot(vv(:,2),'k','linewidth',1.5);title('PC2')
        set(ps(2),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[pyl pyh]);
        ps(3) = subplot(333);
        plot(vv(:,3),'k','linewidth',1.5);title('PC3')
        set(ps(3),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[pyl pyh]);
        % Upper percentile subplots
        ps(4) = subplot(334);
        plot(mean(wsPCahi),'k','linewidth',1.5);title('95th percentile'); hold on
        plot(mean(wsPCahi)+std(wsPCahi),'k','linewidth',1,'linestyle',':');
        plot(mean(wsPCahi)-std(wsPCahi),'k','linewidth',1,'linestyle',':'); hold off
        set(ps(4),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        ps(5) = subplot(335);
        plot(mean(wsPCbhi),'k','linewidth',1.5); title('95th percentile'); hold on
        plot(mean(wsPCbhi)+std(wsPCbhi),'k','linewidth',1,'linestyle',':');
        plot(mean(wsPCbhi)-std(wsPCbhi),'k','linewidth',1,'linestyle',':');hold off
        set(ps(5),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        ps(6) = subplot(336);
        plot(mean(wsPCchi),'k','linewidth',1.5); title('95th percentile'); hold on
        plot(mean(wsPCchi)+std(wsPCchi),'k','linewidth',1,'linestyle',':');
        plot(mean(wsPCchi)-std(wsPCchi),'k','linewidth',1,'linestyle',':'); hold off
        set(ps(6),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        % Lower percentile subplots
        ps(7) = subplot(337);
        plot(mean(wsPCalo),'k','linewidth',1.5); title('5th percentile'); hold on
        plot(mean(wsPCalo)+std(wsPCalo),'k','linewidth',1,'linestyle',':');
        plot(mean(wsPCalo)-std(wsPCalo),'k','linewidth',1,'linestyle',':'); hold off
        set(ps(7),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        ps(8) = subplot(338);
        plot(mean(wsPCblo),'k','linewidth',1.5); title('5th percentile'); hold on
        plot(mean(wsPCblo)+std(wsPCblo),'k','linewidth',1,'linestyle',':');
        plot(mean(wsPCblo)-std(wsPCblo),'k','linewidth',1,'linestyle',':'); hold off
        set(ps(8),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        ps(9) = subplot(339);
        plot(mean(wsPCclo),'k','linewidth',1.5); title('5th percentile'); hold on
        plot(mean(wsPCclo)+std(wsPCclo),'k','linewidth',1,'linestyle',':');
        plot(mean(wsPCclo)-std(wsPCclo),'k','linewidth',1,'linestyle',':'); hold off
        set(ps(9),'color','none','box','on','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
        try waitforbuttonpress; close(pfig); catch; end %  Added close upon key press functionality 
    end

% ORIGIN3R
    function origin3r(source,~)
        set(sb(7), 'Enable', 'off'); drawnow; set(sb(7), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        oval = source.Value;
        if oval == get(source,'Max')
            origins = 0;
            [T1,T2]=view(ms);plott3r; % Update plot
        elseif oval == get(source,'Min')
            origins = 1;
            [T1,T2]=view(ms);plott3r; % Update plot
        end
    end

% COLOR3R
    function color3r(source,~)
        %  Added color3r function to enable changing plot background color.
        set(sb(27), 'Enable', 'off'); drawnow; set(sb(27), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        oval = source.Value;
        if oval == get(source,'Max')
            plotcolor = [1 1 1];
            set(sb(27),'String','Gray background')
        elseif oval == get(source,'Min')
            plotcolor = 'none';
            set(sb(27),'String','White background')
        end
        [T1,T2]=view(ms);plott3r;plott3r2D; % Update plots
    end

% ADD3R
    function add3r(source,~)
        set(sb(8), 'Enable', 'off'); drawnow; set(sb(8), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        if length(pidx)>1;errordlg('Plot one cluster at a time to add additional cluster');end %  Updated cluster/unit terminology
        if isempty(find(cellfun(@isempty,cidx)));errordlg('No vacancies available for additional cluster');end %  Updated cluster/unit terminology
        rotate3d off; pan off; zoom off; datacursormode off; brush off; plotedit off
        
        X0=[];Y0=[];Z0=[];
        set(ms,'UserData',0)
        set(msp,'HitTest','off');
        set(ms,'ButtonDownFcn',@dragg3r0);
        set(gf,'WindowButtonUpFcn',@dragg3r1);
        waitfor(ms,'UserData',1);
        [T1,T2]=view(ms);plott3r;plott3r2D;wav3r;hist3r; % Update plots
        
        function dragg3r0(varargin)
            set(gf,'WindowButtonMotionFcn',@dragg3r);
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            X0 = x0; Y0 = y0; Z0 = z0;
        end
        function dragg3r(varargin)
            hold on
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            X0 = [X0 x0]; Y0 = [Y0 y0]; Z0 = [Z0 z0];
            plot3(X0,Y0,Z0,'k','LineWidth',2,'ButtonDownFcn',@dragg3r0);
            drawnow
        end
        function dragg3r1(varargin)
            set(gf,'WindowButtonMotionFcn','');  % Eliminate fcn on release
            set(gf,'WindowButtonUpFcn','');
            set(ms,'ButtonDownFcn',''); % Limit to 1 ROI
            set(msp,'HitTest','on');
            % Update plot with selected waveforms
            T=view; % Transformation matrix for current view
            % Convert X,Y,Z data coordinates to screen projection
            XYZ = T*[X,Y,Z,ones(nspk,1)]';
            % Convert polygon coordinates to screen projection
            msdata = get(ms,'Children'); % Contour points
            x = get(msdata,'XData'); y = get(msdata,'YData'); z = get(msdata,'ZData');
            xyz = T*[x{1}; y{1}; z{1}; ones(1,length(x{1}))];
            % Assign points inside polygon to next empty cluster, update current cluster
            nidx = min(find(cellfun(@isempty,cidx)));
            awavs = intersect(find(inpolygon(XYZ(1,:),XYZ(2,:),xyz(1,:),xyz(2,:))),cidx{pidx});
            cidx{nidx} = awavs;
            cidx{pidx} = setdiff(cidx{pidx},awavs);
            pidx = nidx;
            set(ms,'UserData',1)
        end
    end

% REMOV3R
    function remov3r(source,~)
        set(sb(9), 'Enable', 'off'); drawnow; set(sb(9), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        if length(pidx)>1;errordlg('Plot one cluster at a time to remove waveforms');end
        rotate3d off; pan off; zoom off; datacursormode off; brush off; plotedit off
        
        X0=[];Y0=[];Z0=[];
        set(ms,'UserData',0)
        set(msp,'HitTest','off');
        set(ms,'ButtonDownFcn',@dragg3r0)
        set(gf,'WindowButtonUpFcn',@dragg3r1);
        waitfor(ms,'UserData',1);
        [T1,T2]=view(ms);
        plott3r;plott3r2D;wav3r;hist3r; % Update plots
        
        function dragg3r0(varargin)
            set(gf,'WindowButtonMotionFcn',@dragg3r);
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            X0 = x0; Y0 = y0; Z0 = z0;
        end
        function dragg3r(varargin)
            hold on
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            X0 = [X0 x0]; Y0 = [Y0 y0]; Z0 = [Z0 z0];
            plot3(X0,Y0,Z0,'k','LineWidth',2,'ButtonDownFcn',@dragg3r0);
            drawnow
        end
        function dragg3r1(varargin)
            set(gf,'WindowButtonMotionFcn','');  % Eliminate fcn on release
            set(gf,'WindowButtonUpFcn','');
            set(ms,'ButtonDownFcn',''); % Limit to 1 ROI
            set(msp,'HitTest','on');
            % Update plot with selected waveforms
            T=view; % Transformation matrix for current view
            % Convert X,Y,Z data coordinates to screen projection
            XYZ = T*[X,Y,Z,ones(nspk,1)]';
            % Convert contour coordinates to screen projection
            msdata = get(ms,'Children'); % Contour points
            x = get(msdata,'XData'); y = get(msdata,'YData'); z = get(msdata,'ZData');
            xyz = T*[x{1}; y{1}; z{1}; ones(1,length(x{1}))];
            % Assign points inside polygon to Unsorted, update current cluster
            rwavs = intersect(find(inpolygon(XYZ(1,:),XYZ(2,:),xyz(1,:),xyz(2,:))),cidx{pidx});
            cidx{1} = sort([cidx{1} rwavs]);
            cidx{pidx} = setdiff(cidx{pidx},rwavs);
            set(ms,'UserData',1)
        end
    end

% EXPLOR3R
    function explor3r(source,~)
        set(sb(10), 'Enable', 'off'); drawnow; set(sb(10), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        rotate3d off; pan off; zoom off; datacursormode off; brush off; plotedit off
        esmpl = 50; persistent ef
        set(ms,'ButtonDownFcn',@click3r)
        set(msp,'HitTest','off');
        
        function click3r(varargin)
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            T=view; % Transformation matrix for current view
            XYZ = T*[X,Y,Z,ones(nspk,1)]';
            xyz = T*[x0; y0; z0; 1];
            edvec = (XYZ(1,:)-xyz(1,:)).^2+(XYZ(2,:)-xyz(2,:)).^2;
            edvec0=sort(edvec); edvals=edvec0(1:esmpl);
            % Label points
            if length(ms.Children)>esmpl;delete(ms.Children(1:esmpl)); end
            hold on
            for eii = 1:esmpl
                wid = find(edvec==edvals(eii));
                plot3(X(wid),Y(wid),Z(wid),'ro','LineWidth',1.5);
            end
            hold off
            % Plot waves
            if ishandle(ef); delete(ef);end
            ef = figure('units','normalized','position',[.6,.4,.3,.4],'toolbar','none','menu','none','name','','numbertitle','off');
            title('[CLOSE WINDOW TO EXIT]'); %  added instruction title to clarify user must close explorer subplot window to exit function 
            hold on
            for eii = 1:esmpl
                wid = find(edvec==edvals(eii));
                plot(waves(wid,:),'linewidth',1.5,'color',colorvec(1,:));
            end
            hold off
            set(ef.CurrentAxes,'color','none','box','off','xcolor','none','ycolor','none','xtick',[],'ytick',[],'xlim',[1 wxh],'ylim',[wyl wyh]);
            set(ef,'CloseRequestFcn',@click3r1);
        end
        function click3r1(varargin)
            set(ms,'ButtonDownFcn',''); % Limit to 1 ROI
            set(msp,'HitTest','on');
            delete(ef);
            delete(ms.Children(1:esmpl));
        end
    end

% EXCIS3R
    function excis3r(source,~)
        set(sb(11), 'Enable', 'off'); drawnow; set(sb(11), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        if length(pidx)>1 || ~isempty(find(~cellfun(@isempty,cidx(2:end)))) ;errordlg('Excise waveforms can only be performed when all waves are assigned to Cluster 1');end %  Updated to avoid crash when excis3r is used with more than one cluster defined 
        rotate3d off; pan off; zoom off; datacursormode off; brush off; plotedit off
        
        X0=[];Y0=[];Z0=[];
        set(ms,'UserData',0)
        set(msp,'HitTest','off');
        set(ms,'ButtonDownFcn',@dragg3r0)
        set(gf,'WindowButtonUpFcn',@dragg3r1);
        waitfor(ms,'UserData',1);
        
        [T1,T2]=view(ms);featur3r;plott3r;plott3r2D;wav3r;hist3r; % Update plots
        
        function dragg3r0(varargin)
            set(gf,'WindowButtonMotionFcn',@dragg3r);
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            X0 = x0; Y0 = y0; Z0 = z0;
        end
        function dragg3r(varargin)
            hold on
            pt = get(ms,'CurrentPoint');
            x0 = pt(1,1); y0 = pt(1,2); z0 = pt(1,3);
            X0 = [X0 x0]; Y0 = [Y0 y0]; Z0 = [Z0 z0];
            plot3(X0,Y0,Z0,'k','LineWidth',2,'ButtonDownFcn',@dragg3r0);
            drawnow
        end
        function dragg3r1(varargin)
            set(gf,'WindowButtonMotionFcn','');  % Eliminate fcn on release
            set(gf,'WindowButtonUpFcn','');
            set(ms,'ButtonDownFcn',''); % Limit to 1 ROI
            set(msp,'HitTest','on');
            % Update plot with selected waveforms
            T=view; % Transformation matrix for current view
            % Convert X,Y,Z data coordinates to screen projection
            XYZ = T*[X,Y,Z,ones(nspk,1)]';
            % Convert contour coordinates to screen projection
            msdata = get(ms,'Children'); % Contour points
            x = get(msdata,'XData'); y = get(msdata,'YData'); z = get(msdata,'ZData');
            xyz = T*[x{1}; y{1}; z{1}; ones(1,length(x{1}))];
            % Remove points inside polygon from waves and spiketimes
            ewavs = intersect(find(inpolygon(XYZ(1,:),XYZ(2,:),xyz(1,:),xyz(2,:))),cidx{pidx});
            waves = waves(setdiff(cidx{pidx},ewavs),:);
            spiketimes = spiketimes(setdiff(cidx{pidx},ewavs));
            set(ms,'UserData',1)
        end
    end

% ALLMATCH3R
    function allmatch3r(source,~)
        set(sb(16), 'Enable', 'off'); drawnow; set(sb(16), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        tidx = find(~cellfun(@isempty,cidx));
        if length(tidx)>1
            % cluster templates
            for tii=1:length(tidx)
                templates(tii,:) = mean(waves(cidx{tidx(tii)},:));
            end
            % Euclidean distances
            for tjj = 1:length(waves)
                for tii = 1:length(tidx)
                    eucdists(tjj,tii) = sqrt(sum((waves(tjj,:)-templates(tii,:)).^2));
                end
            end
            [~,eucvec] = min(eucdists,[],2);
            % Update cidx
            for tii=1:length(tidx)
                cidx{tidx(tii)} = find(eucvec==tii)';
            end
            pidx = tidx; % Plot all clusters included in template matching
            [T1,T2]=view(ms);plott3r;plott3r2D;wav3r;hist3r; % Update plots
        end
    end

% COMBIN3R
    function combin3r(source,~)
        set(sb(12), 'Enable', 'off'); drawnow; set(sb(12), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        cfw=200; cfh=200;
        cfig = figure('units','pixels','position',[100,100,cfw,cfh],'toolbar','none','menu','none');
        set(cfig,'name','Combine which?','numbertitle','off');
        % cluster checkboxes
        ccb(1) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.86,70,15],'string','Cluster 1');
        ccb(2) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.70,70,15],'string','Cluster 2');
        ccb(3) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.54,70,15],'string','Cluster 3');
        ccb(4) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.38,70,15],'string','Cluster 4');
        ccb(5) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.22,70,15],'string','Cluster 5');
        ccb(6) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.06,70,15],'string','Cluster 6');
        % Combine button
        uicontrol('style','pushbutton','units','pixels',...
            'position',[100,cfh*.45,70,20],'string','Combine',...
            'callback',@check3r);
        function check3r(varargin)
            ccbs = get(ccb,'Value');
            checked = find([ccbs{:}]);
            if length(checked) > 1
                for cii = 2:length(checked)
                    cidx{checked(1)} = sort([cidx{checked(1)} cidx{checked(cii)}]);
                    cidx{checked(cii)}=[];
                end
                for cii = 1:length(checked)
                    set(ub4(checked(cii)),'String',''); %  Reset unit number strings after combining
                end
                close(cfig);
                pidx = min(checked);
                [T1,T2]=view(ms);plott3r;plott3r2D;wav3r;hist3r; % Update plots
            else
                close(cfig);
            end
        end
    end

% CHANNEL3R
    function channel3r(source,~)
        set(sb(20), 'Enable', 'off'); drawnow; set(sb(20), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        cfw=200; cfh=200;
        cfig = figure('units','pixels','position',[150,100,cfw,cfh],'toolbar','none','menu','none');
        set(cfig,'name','Include which channels?','numbertitle','off');
        % cluster checkboxes
        ccb(1) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.86,70,15],'string','Channel 1');
        ccb(2) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.70,70,15],'string','Channel 2');
        ccb(3) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.54,70,15],'string','Channel 3');
        ccb(4) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.38,70,15],'string','Channel 4');
        ccb(5) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,cfh*.22,70,15],'string','Channel 5');
        % Combine button
        uicontrol('style','pushbutton','units','pixels',...
            'position',[100,cfh*.45,70,20],'string','Enter',...
            'callback',@check3r);
        function check3r(varargin)
            ccbs = get(ccb,'Value');
            checked = find([ccbs{:}]);
            if ~isempty(checked)
                nWavSamples = opts.nWavSamples;
                keepWavs = [];
                for cii = 1:length(checked)
                    keepWavs = [keepWavs 1+(checked(cii)*nWavSamples)-nWavSamples:checked(cii)*nWavSamples];
                end
                waves = waves(:,keepWavs);
                close(cfig);
                [T1,T2]=view(ms);featur3r;plott3r;plott3r2D;wav3r;hist3r; % Update plots
            else
                close(cfig);
            end
        end
    end

% PIDX3R
    function pidx3r(source,~)
        set(sb(13), 'Enable', 'off'); drawnow; set(sb(13), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        pfw=200; pfh=200;
        pfig = figure('units','pixels','position',[100,100,pfw,pfh],'toolbar','none','menu','none');
        set(pfig,'name','Plot which?','numbertitle','off');
        % cluster checkboxes
        pcb(1) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,pfh*.86,70,15],'string','Cluster 1');
        pcb(2) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,pfh*.70,70,15],'string','Cluster 2');
        pcb(3) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,pfh*.54,70,15],'string','Cluster 3');
        pcb(4) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,pfh*.38,70,15],'string','Cluster 4');
        pcb(5) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,pfh*.22,70,15],'string','Cluster 5');
        pcb(6) = uicontrol('style','checkbox','units','pixels',...
            'position',[10,pfh*.06,70,15],'string','Cluster 6');
        % Plot button
        uicontrol('style','pushbutton','units','pixels',...
            'position',[100,pfh*.5,70,20],'string','Plot',...
            'callback',@check3r);
        function check3r(varargin)
            pcbs = get(pcb,'Value');
            checked = find([pcbs{:}]);
            if ~isempty(checked)
                pidx = checked;
                close(pfig);
                [T1,T2]=view(ms);plott3r; % Update plots
            else
                close(pfig);
            end
        end
    end

% ALLPIDX3R
    function allpidx3r(source,~)
        set(sb(14), 'Enable', 'off'); drawnow; set(sb(14), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        pidx = find(~cellfun(@isempty,cidx)); %  Updated to fix add3r bug when only one cluster is present
        [T1,T2]=view(ms);plott3r; % Update plots
    end

% XCORR3R
    function xcorr3r(source,~)
        plott3r;plott3r2D;wav3r;hist3r;
        set(sb(17), 'Enable', 'off'); drawnow; set(sb(17), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        xfig = figure('units','normalized','outerposition',[.25,.25,.5,.5],'toolbar','none','menu','none');
        set(xfig,'name','Auto- and cross-correlation matrices','numbertitle','off');
        xidx = find(~cellfun(@isempty,cidx)); xsz = length(xidx);
        for xii = 1:length(xidx)
            for xjj = 1:length(xidx)
                tres = 1; % ms  %AH Increased
                maxlag = 40; % ms  %ah Increased - too low to ID burstiness
                switch opts.tspkUnits
                    case 's'
                        spk1 = spiketimes(cidx{xii})*1000;
                        spk2 = spiketimes(cidx{xjj})*1000;
                    case 'ms'
                        spk1 = spiketimes(cidx{xii});
                        spk2 = spiketimes(cidx{xjj});
                    case 'samples'
                        spk1 = (spiketimes(cidx{xii})/opts.Fs)*1000;
                        spk2 = (spiketimes(cidx{xjj})/opts.Fs)*1000;
                end
                binvec = 0:tres:max([spk1; spk2]);
                raster1 = histc(spk1,binvec)./tres;
                raster2 = histc(spk2,binvec)./tres;
                [c,lags] = xcorr(raster1,raster2,maxlag);
                subplot(xsz,xsz,xii*xsz-(xsz-xjj))
                if xii==xjj
                    xcolor = colorvec(xii,:);
                    c(find(lags==0)) = NaN;
                else
                    xcolor = [.5,.5,.5];
                end
                bar(lags.*tres,c,'FaceColor',xcolor,'EdgeColor','none','BarWidth',1);
                % set(gca,'color','none','xcolor','none','ycolor','none','xtick',[],'ytick',[],'xlim',[-.5-maxlag/2 +.5+maxlag/2],'ylim',[min(c)-(.1*(max(c)-min(c))) max(c)+(.1*(max(c)-min(c)))]);
                set(gca,'box','off','color','none','ycolor','none','TickDir','out','xtick',-maxlag:5:maxlag,'ytick',[],'xlim',[-.5-maxlag/2 +.5+maxlag/2],'ylim',[0 max(c)+(.1*(max(c)-min(c)))],'fontsize',fntsz);
                drawnow
            end
        end
        try waitforbuttonpress; close(xfig); catch; end %  Added close upon key press functionality 
    end

% RESETT3R
    function resett3r(source,~)
        set(sb(22), 'Enable', 'off'); drawnow; set(sb(22), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        waves = waves0; spiketimes = spiketimes0;
        cidx{1} = 1:length(spiketimes);
        for rii = 2:length(cidx)
            cidx{rii}=[];
            set(ub4(rii),'String',''); %  Reset unit ID string 
        end
        pidx = 1;
        [T1,T2]=view(ms);featur3r;plott3r;plott3r2D;wav3r;hist3r; % Update plots
    end

% SOLOPIDX3R
    function solopidx3r(source,~,keyID) % Added third input arg for calling with hotk3ys
        if nargin < 3
            pidx = source.UserData;
        else
            pidx = keyID;
        end
        [T1,T2]=view(ms);plott3r; % Update plots
        set(ub1(pidx), 'Enable', 'off'); drawnow; set(ub1(pidx), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
    end

% DELET3R
    function delet3r(source,~)
        uid = source.UserData;
        if uid~=1
            cidx{1} = sort([cidx{1} cidx{uid}]);
            cidx{uid}=[]; pidx = 1;
            [T1,T2]=view(ms);plott3r;plott3r2D;wav3r;hist3r; % Update plots
        end
        set(ub2(uid), 'Enable', 'off'); drawnow; set(ub2(uid), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        set(ub4(uid),'String',''); %  Reset unit ID string 
    end

% ALLDELET3R
    function alldelet3r(source,~)
        set(sb(15), 'Enable', 'off'); drawnow; set(sb(15), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        cidx{1} = 1:length(spiketimes);
        for rii = 2:length(cidx)
            cidx{rii}=[];
            set(ub4(rii),'String',''); %  Reset unit ID string
        end
        pidx = 1;
        [T1,T2]=view(ms);plott3r;plott3r2D;wav3r;hist3r; % Update plots %  removed featur3r call to avoid unnecessary recalculating features
    end

% MATCH3R
    function match3r(source,~)
        tid = source.UserData;
        tidx = find(~cellfun(@isempty,cidx));
        tidx = tidx(tidx~=tid);
        if ~isempty(tidx)
            % cluster templates
            for tii=1:length(tidx)
                templates(tii,:) = mean(waves(cidx{tidx(tii)},:));
            end
            % Euclidean distances
            twaves = waves(cidx{tid},:);
            for tjj = 1:length(twaves)
                for tii = 1:length(tidx)
                    eucdists(tjj,tii) = sqrt(sum((twaves(tjj,:)-templates(tii,:)).^2));
                end
            end
            [~,eucvec] = min(eucdists,[],2);
            % Update cidx
            for tii=1:length(tidx)
                cidx{tidx(tii)} = sort([cidx{tidx(tii)} cidx{tid}(eucvec==tii)]);
            end
            cidx{tid}=[];
            pidx = tidx; % Plot all clusters included in template matching
            [T1,T2]=view(ms);plott3r;plott3r2D;wav3r;hist3r; % Update plots
        end
        set(ub3(tid), 'Enable', 'off'); drawnow; set(ub3(tid), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
    end

% VIEW3R
    function view3r(source,~,keyID) %  Added third input arg for calling with hotk3ys
        if nargin < 3
            vid = source.UserData;
        else
            vid = keyID;
        end
        set(tb(vid), 'Enable', 'off'); drawnow; set(tb(vid), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        switch vid
            case 1
                Xid = 'Spiketimes';     X = eval(Xid);
                Yid = 'PC1';            Y = eval(Yid);
                Zid = 'PC2';            Z = eval(Zid);
                T1=90;T2=0;
                plott3r; % Update plots
            case 2
                Xid = 'Spiketimes';     X = eval(Xid);
                Yid = 'PC1';            Y = eval(Yid);
                Zid = 'PeakValley';     Z = eval(Zid);
                T1=90;T2=0;
                plott3r; % Update plots
            case 3
                Xid = 'Spiketimes';     X = eval(Xid);
                Yid = 'PC2';            Y = eval(Yid);
                Zid = 'PC1';            Z = eval(Zid);
                T1=0;T2=0;
                plott3r; % Update plots
            case 4
                Xid = 'Spiketimes';     X = eval(Xid);
                Yid = 'PC1';            Y = eval(Yid);
                Zid = 'PeakValley';     Z = eval(Zid);
                T1=0;T2=0;
                plott3r; % Update plots
        end
        % Update popup menu strings (JB updated 2018-08-13)
        fidx2 = find( strcmp( fstrings, Xid) ); set( sb(2), 'Value', fidx2);
        fidx4 = find( strcmp( fstrings, Yid) ); set( sb(4), 'Value', fidx4);
        fidx6 = find( strcmp( fstrings, Zid) ); set( sb(6), 'Value', fidx6);
    end

% SAV3R
    function sav3r(source,~)
        
        if max(diff(find(~cellfun(@isempty,cidx))))>1;errordlg('Set clusters to lowest available cluster slot'); end % Added to avoid crash in cluster_stats, e.g., if cluster slots 1, 2, and 4 are populated but 3 is empty 
        
        %  Close any figures other than the GUI (e.g., unclosed utility subplots)
        all_figs = findobj(0, 'type', 'figure');
        for fg = 1:numel(all_figs)
            if ~strcmp( all_figs(fg).Tag, 'GUI' )
                set(all_figs(fg),'CloseRequestFcn','')
                delete( all_figs(fg) );
            end
        end
        
        %  Make sure unit sort quality codes are entered
        missingValues = 0;
        for uz = 1:length(cidx)
            if ~isempty(cidx{uz}) && strcmp('', get(ub4(uz),'String'))
                missingValues = 1;
                break
            end
        end
        if sum(missingValues)>0
            % warning('Sort codes missing'); %  Updated for consisency with GUI terminology
        end
        
        % Save figure and data files unless existing sort was loaded and no edits were made
        if existing_sort_flag && isequal(cidx, cidx_existing)
            saveflag = false;
            set(gf,'UserData',1)            
        else
            saveflag = true;
            set(gf,'UserData',1)
        end
        
    end

% NORM3R
    function [X] = norm3r(X0)
        X = (X0-min(X0))/(max(X0)-min(X0));
    end

% SPKISI3R
    function [bins,isi1] = spkisi3r(tspk)
        switch opts.tspkUnits
            case 's'
                isi = diff(tspk)*1000;
            case 'ms'
                isi = diff(tspk);
            case 'samples'
                isi = diff(tspk/opts.Fs)*1000;
        end
        isi1 = 100 * length(find(isi<3))./length(isi); % Short ISI count
        isiBins=logspace(-1,4,60);
        bins=histc(isi,isiBins);
    end

% FIGUR3R
    function figur3r(source,~, final)
        set(sb(21), 'Enable', 'off'); drawnow; set(sb(21), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        allpidx3r;
        
        % Plot stats
        %{
        if length(~cellfun('isempty',cidx)) > 1
            [~,fhan] = clust3r_stats(cidx(2:max(find(~cellfun(@isempty,cidx)))),[Y Z],waves,spiketimes,1,colorvec);
            drawnow;
            % pause(3); % may take a few secs for cluster stats to finish calculating and plotting
        end
        %}
        
        fign = fign + 1;
        % Give the final figure a zero
        %         if fign > 1
        if exist('final', 'var') && final
            figstring = [figpath filesep 'Channel_' num2str(chan,'%02d') '_final.png']; % designate 'final' figure when user selects 'Save & Exit'
        else
            figstring = [figpath filesep 'Channel_' num2str(chan,'%02d') '_supplemental_' num2str(fign) '.png']; % designate 'supplemental' figure when user selects 'Save figure'
        end
        
        imageData = frame2im(getframe(gf));
        if exist('fhan', 'var') && isa(fhan, 'matlab.ui.Figure')
            imageData0 = frame2im(getframe(fhan));
            h_diff = size(imageData,1) - size(imageData0,1);
            imageData0 = [ zeros(h_diff, size(imageData0,2), size(imageData0,3))+240 ; imageData0];
            imageData = [imageData imageData0];
            close(fhan);
        end
        imwrite(imageData,figstring)

    end

% STAT3R
    function stat3r(source,~)
        plott3r;plott3r2D;wav3r;hist3r;
        set(sb(18), 'Enable', 'off'); drawnow; set(sb(18), 'Enable', 'on'); % Return focus to GUI for enabling hotkeys
        if length(find(~cellfun(@isempty,cidx)))<2;errordlg('Cluster statistics require at least two clusters');end
        [cluster_stats_n,fhan,savestats] = clust3r_stats(cidx(2:max(find(~cellfun(@isempty,cidx)))),[Y Z],waves,spiketimes,1,[],colorvec); % added fhan output for close upon key press functionality 
        if ~isempty(cluster_stats_n)
            cluster_stats_n.features{1}= Yid;
            cluster_stats_n.features{2}= Zid;
            for sii = 1:length(cidx)
                if ~isempty(cidx{sii})
                    cluster_stats_n.qual{sii} = get(ub4(sii),'String');
                end
            end
            if savestats
                if isempty(cluster_statistics)
                    cluster_statistics = cluster_stats_n;
                else
                    cluster_statistics(end+1) = cluster_stats_n;
                end
            end
        end
        try waitforbuttonpress; close(fhan); catch; end %  Added close upon key press functionality
    end

% HOTK3YS
    function hotk3ys(source, e) %  added function for hotkey mappings
        switch e.Key
            
                %  Side buttons  
            case 'f'
                figur3r;
            case 'd'
                alldelet3r;
            case 'e'
                excis3r;
            case 'p'
                allpidx3r;
            case 'a'
                add3r;
            case 'r'
                remov3r;
            case 'c'
                stat3r;
            case 'x'
                xcorr3r;
            case 't'
                allmatch3r;
            case 'return'
                sav3r;
             case 'q' % Quit without saving
                exit3r
                
                %  Enable cluster views by keyboard or numpad  
            case '1'
                solopidx3r([],[],1)
            case '2'
                solopidx3r([],[],2)
            case '3'
                solopidx3r([],[],3)
            case '4'
                solopidx3r([],[],4)
            case '5'
                solopidx3r([],[],5)
            case '6'
                solopidx3r([],[],6)
            case 'numpad1'
                solopidx3r([],[],1)
            case 'numpad2'
                solopidx3r([],[],2)
            case 'numpad3'
                solopidx3r([],[],3)
            case 'numpad4'
                solopidx3r([],[],4)
            case 'numpad5'
                solopidx3r([],[],5)
            case 'numpad6'
                solopidx3r([],[],6)
                
                %  2D plot views
            case 'f1'
                view3r([],[],1)
            case 'f2'
                view3r([],[],2)
            case 'f3'
                view3r([],[],3)
            case 'f4'
                view3r([],[],4)
                
        end
    end

% EXIT3R
    function exit3r(source,~) %  Added exit3r function for exiting without saving
        
        opts.Interpreter = 'tex';
        opts.Default = 'Yes';
        quest = 'Exit without saving current channel?';
        answer = questdlg(quest,'Confirm exit',...
            'Yes','No',opts);
        
        if strcmp( answer, 'Yes')
            %  Close any figures other than the GUI (e.g., unclosed utility subplots)
            all_figs = findobj(0, 'type', 'figure');
            for fg = 1:numel(all_figs)
                if ~strcmp( all_figs(fg).Tag, 'GUI' )
                    set(all_figs(fg),'CloseRequestFcn','')
                    delete( all_figs(fg) );
                end
            end
            rstrings{end+1} = 'return'; 
            saveflag = false;
            set(gf,'UserData',1)
        else
        end

    end

end


