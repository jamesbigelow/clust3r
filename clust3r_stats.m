function [stats, fig_han, savestats ] = clust3r_stats( cidx, wave_features, waveforms, spiketimes, plotflag, gui_flag, colors)

% NAME: CLUST3R_STATS
%
% DESCRIPTION: Aux function for clust3r spike-sorting GUI which calculates
% cluster isolation statistics with optional plot. Incorporates several 
% utilities adapted from KFMMAutosorter written by Matthew Fellows.
%
% INPUT:
% cidx: Cluster indices indicating which spikes belong to which cluster, N length vector 
% wave_features: Spike waveform features used for spike sorting, 2 x N array (two features x n spikes)
% waveforms: Spike waveforms, N x K array (spike index x n samples)
% spiketimes: Spiketimes, N length vector (default sec; see params)
% plotflag: True/False generate plot
% gui_flag: True/False generate plot as GUI for optional save stats function
% colors: cluster colors 
% 
% OUTPUT:
% stats: Structure containing cluster isolation statistics 
% fig_han: Figure handle
% savestats: True/False store stats in final clust3r output files (set by GUI button)
% 
% INTERNAL FUNCTIONS: 
% fun1, get_l_ratio, get_isolation_distance, 

fntsz = 10; if max(get(groot,'Screensize'))<2000; fntsz = fntsz*.6; end
if ~exist('gui_flag','var')
    gui_flag = false;
end

cidx_tmp = ones(1,length(spiketimes));
for ii = 1:length(cidx)
    cidx_tmp( cidx{ii}) = ii+1;
end

silcenters = -.5:.05:1;
r = randperm(numel(cidx_tmp),min(500,numel(cidx_tmp)));
cidx_unique = unique(cidx_tmp(r));
stats = struct();
fig_han = [];

if numel(cidx_unique) > 1
    
    uids = unique(cidx_tmp);
    
    %% Short ISI
    for ux = 1:length(uids)
        tspk = spiketimes(find(cidx_tmp == uids(ux)));
        isi = diff(tspk);
        short_isi_1_ms(ux) = 100 * length(find(isi<1))./length(isi);
        short_isi_3_ms(ux) = 100 * length(find(isi<3))./length(isi);
    end
    stats.short_isi_1_ms = short_isi_1_ms;
    stats.short_isi_3_ms = short_isi_3_ms;
    
    %% Waveform SNR
    wv_mu = waveforms(find(cidx_tmp == 1),:);
    for ux = 1:length(uids)
        wv_su = waveforms(find(cidx_tmp == uids(ux)),:);
        wave_snr(ux) =  ( max(mean(wv_su)) - min(mean(wv_su)) ) / ( max(mean(wv_mu)) - min(mean(wv_mu)) );
    end
    stats.wave_snr = wave_snr;
    
    %% Silhouette indices
    % PCA based ---------
    s = silhouette(wave_features(r,:),cidx_tmp(r));
    sils.Feature.avg = mean(s);
    for ii = 1:numel(cidx_unique)
        ss(ii) = mean(s(find(cidx_tmp(r) == cidx_unique(ii))));
        if(plotflag)
            hh = hist(s(find(cidx_tmp(r) == cidx_unique(ii))),silcenters);
            hh = hh/sum(hh);
            hhs(ii,:) = hh;
            labels{ii} = sprintf('C%d',cidx_unique(ii));
            sils.Feature.hhs = hhs;
        end
    end
    sils.Feature.clustavgs = ss;
    sils.Feature.cidx_unique = cidx_unique;
    sils.Feature.cidx_tmp = cidx_tmp;
    sils.Feature.s = s;
    
    % Wave based ---------
    t = silhouette(waveforms(r,:),cidx_tmp(r));
    sils.Wave.avg = mean(t);
    for ii = 1:numel(cidx_unique)
        ss(ii) = mean(t(find(cidx_tmp(r) == cidx_unique(ii))));
        if(plotflag)
            h = hist(t(find(cidx_tmp(r) == cidx_unique(ii))),silcenters);
            hh = hh/sum(hh);
            hhs(ii,:) = hh;
            labels{ii} = sprintf('C%d',cidx_unique(ii));
            sils.Wave.hhs = hhs;
        end
    end
    sils.Wave.clustavgs = ss;
    sils.Wave.inds = cidx_unique;
    sils.Wave.clusterBy = cidx_tmp;
    sils.Wave.s = t;
    
    %% Isolation distance and L Ratio
    for ii = 1:numel(cidx_unique)
        cidx_tmp0 = find(cidx_tmp == ii);
        [IsolDist(ii), mNoise{ii}, mCluster{ii}] = get_isolation_distance(wave_features,cidx_tmp0);
        [P{ii}, L(ii), Lratio(ii), df(ii), n_spk_clust0(ii)] = get_l_ratio(wave_features, cidx_tmp0);
    end
    
    idist.isoldist = IsolDist;
    idist.mNoise = mNoise;
    idist.mCluster = mCluster;
    idist.cidx_tmp = cidx_tmp;
    lrat.LRatio = Lratio;
    lrat.P = P;
    lrat.L = L;
    lrat.df = df;
    lrat.nclusterspikes = n_spk_clust0;
    lrat.cidx_tmp = cidx_tmp;
    
    stats.sils = sils;
    stats.idist = idist;
    stats.lrat = lrat;
    
    if plotflag
        
        fig_han = figure('units','normalized','outerposition',[.8 .05 .2 .45],'toolbar','none','menu','none');
        set(fig_han,'name','Cluster statistics','numbertitle','off')
        sb = uicontrol('Style', 'togglebutton', 'String', 'Save stats','units','normalized','Position',[.75 .9 .15 .05],'Callback',@fun1);
        savestats = get(sb,'value');
        
        if ~exist('colors','var')
            colors = [.2,.2,.2; 1,0,0; 0,0,1; 1,0,1; 0,1,1; 0,1,0;];
        end
        subplot(321)
        for ii = 1:numel(cidx_unique)
            plot(silcenters,sils.Feature.hhs(ii,:)','color',colors(ii,:),'linewidth',2)
            hold on
            txtdsc{ii} = sprintf('Cluster %d:       %2.2f\n', cidx_unique(ii), sils.Feature.clustavgs(ii) );
        end
        axis tight
        set(gca,'color','none');
        lh = legend(labels);
        set(lh,'box','off','location','best','fontsize',fntsz)
        box off
        set(gca,'fontsize',fntsz);
        title('Silhouette (PCA based)','fontsize',fntsz*1.2)
        
        subplot(322)
        tx1 = sprintf('Mean:      %2.2f\n',sils.Feature.avg);
        txt = [tx1,[txtdsc{:}]];
        axis off
        text(0,.3,txt,'fontsize',fntsz*1.2)
        subplot(323)
        set(gca,'fontsize',fntsz)
        for ii = 1:numel(cidx_unique)
            plot(silcenters,sils.Wave.hhs(ii,:)','color',colors(ii,:),'linewidth',2)
            hold on
            txtdsc{ii} = sprintf('Cluster %d:       %2.2f\n', cidx_unique(ii), sils.Wave.clustavgs(ii) );
        end
        axis tight
        set(gca,'color','none');
        lh = legend(labels);
        set(lh,'box','off','location','best','fontsize',fntsz)
        box off
        set(gca,'fontsize',fntsz);
        title('Silhouette (Waveform based)','fontsize',fntsz*1.2)
        
        subplot(324)
        tx1 = sprintf('Mean:      %2.2f\n',sils.Wave.avg);
        txt = [tx1,[txtdsc{:}]];
        axis off
        text(0,.3,txt,'fontsize',fntsz*1.2)
        subplot(325)
        set(gca,'fontsize',fntsz)
        labelnum=0;
        for ii = 1:numel(cidx_unique)
            if isnan(IsolDist(ii))
                txtdsc{ii} = sprintf('Cluster %d isolation      %2.2f\nL Ratio      %2.2f\n',cidx_unique(ii),IsolDist(ii),Lratio(ii));
            else
                labelnum = labelnum + 1;
                labels{2*labelnum-1}=sprintf('Cluster%d',cidx_unique(ii));
                labels{2*labelnum}=sprintf('Noise%d',cidx_unique(ii));
                minyval = min(max(mCluster{ii}),max(mNoise{ii}));
                toplotClust = sort(mCluster{ii}(mCluster{ii}<=minyval));
                toplotNoise = sort(mNoise{ii}(mNoise{ii}<=minyval));
                semilogx(toplotClust,'color',colors(ii,:),'linewidth',2); hold on;
                semilogx(toplotNoise,':','color',colors(ii,:),'linewidth',1); hold on;
                txtdsc{ii} = sprintf('Cluster %d isolation      %2.2f\nL Ratio      %2.2f\n',cidx_unique(ii),IsolDist(ii),Lratio(ii));
            end
        end
        axis tight
        set(gca,'color','none');
        lh = legend(labels);
        ylabel('Mahalanobis distance','fontsize',fntsz)
        xlabel('Cumulative points')
        set(lh,'box','off','location','best','fontsize',fntsz)
        box off
        set(gca,'fontsize',fntsz);
        title('Isolation Distance','fontsize',fntsz,'fontsize',fntsz*1.2)
        subplot(326)
        set(gca,'fontsize',fntsz)
        axis off
        text(0,.3,txtdsc,'fontsize',fntsz*1.2)
        
        if gui_flag
            waitfor(fig_han)
        end
    end
    
else
    stats = struct();
    fig_han = [];
    savestats = [];
end

    function fun1(hObject,~)
        savestats = get(hObject,'Value');
    end


    function [ P, L, Lratio, df, n_spk_clust ] = get_l_ratio( Features, spike_cluster_indices )
        
        n_spk = size(Features,1);
        n_spk_clust = length(spike_cluster_indices);
        unsorted_spikes = setdiff(1:n_spk, spike_cluster_indices);
        
        m = mahal(Features, Features(spike_cluster_indices,:));
        
        df = size(Features,2);
        P = (1-chi2cdf(m(unsorted_spikes),df));
        L = sum(P);
        Lratio = L/n_spk_clust;
        
    end


    function [ IsolDist, mNoise, mCluster ] = get_isolation_distance( Features, spike_cluster_indices )
        
        n_spk = size(Features,1);
        n_spk_clust = length(spike_cluster_indices);
        unsorted_spikes = setdiff(1:n_spk, spike_cluster_indices);
        m = mahal(Features, Features(spike_cluster_indices,:));
        mCluster = m(spike_cluster_indices);
        mNoise = m(unsorted_spikes);
        
        if length(spike_cluster_indices) >= size(Features,2) && length(spike_cluster_indices)<=n_spk/2
            S = sort(mNoise);
            IsolDist = S(n_spk_clust);
        else
            IsolDist = NaN;
        end
        
    end

end