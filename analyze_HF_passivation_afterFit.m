clear all; close all; clc; 

%constant for all samples
T = 25+273.15; 
type = 'p';

%histogram number of bins
nbins = 7; 

%excess carrier densities for evaluation
deltans = [3e14 6e14 8e14 1e15 5e15 1e16]; 

dir='C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set a\3761830s\lifetime spectroscopy\finalized defect fits';
dir1 = [dir '\49a (fired 0 min)'];
dir2 = [dir '\53a (fired 10 min)'];
dir3 = [dir '\54a (unfired 10 min)'];
dir4 = [dir '\61a (unfired 0 min)'];
dir5 = [dir '\60a (control)'];

samp_dir = {'49a','53a','54a','61a','60a';dir1,dir2,dir3,dir4,dir5}; 
[rows,num_samp] = size(samp_dir); 

%load the data we'll need for calculation
load('C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set a\3761830s\set-a_all_data_seta_3761830s.mat'); 
load('C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set a\3761830s\lifetime spectroscopy\defect_fits.mat'); 

for i = 1:num_samp
    %Make the Excel filename to write this data
    Etdatafile = [samp_dir{2,i} '\Etk_plotting.xlsx']; 
    %load the excel data for this sample
    data = xlsread([samp_dir{2,i} '\IDLS two and three curve fitting.xlsm'],'summary','A:G'); 
    %Now we need the doping level and temperature for this sample
    index = find(strcmp(samples,samp_dir{1,i})==1); 
    doping = doping_all{index}(1); %take the first one and assume it's correct
    %get some sample parameters that we'll need later
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,doping,type); 
    %figure out where the lifetime data for this sample is stored
    indextau = find(strcmp(lifetime_analysis(1,:),samp_dir{1,i})==1); 
    defect_now = defect_summary{indextau};
    %Now we need to loop over time
    [times,columns] = size(data); 
    [timestau,columns] = size(defect_now); 
    cm = colormap(hsv(times));
    if times ~= timestau
        disp(['the lifetime and fit times are not the same for sample ' samp_dir{1,i}]); 
    end
    Et_now = cell(times,2); k_now = cell(times,2); alphanN_now = cell(times,2); 
    allEk = figure('units','normalized','outerposition',[0 0 1 1]); 
    dominance = zeros(times,length(deltans)); 
    simple_summary = zeros(times,4); 
    for j = 1:times
        fits_now = data(j,2:5); 
        %defect 1
        [Et_now{j,1},k_now{j,1},alphanN_now{j,1}]=generate_Ek(fits_now(1:2),T,doping,type);
        %defect 2
        [Et_now{j,2},k_now{j,2},alphanN_now{j,2}]=generate_Ek(fits_now(3:4),T,doping,type);
        %write the Et, k, and alphanN for this time to an Excel sheet for
        %easy plotting later
        to_write = [Et_now{j,1}',k_now{j,1}',alphanN_now{j,1}',k_now{j,2}',alphanN_now{j,2}']; 
        labels = {'Et-Ei [eV]','defect 1 k [-]','defect 1 alphanN [1/us]','defect 2 k [-]','alphanN [1/us]'}; 
        xlswrite(Etdatafile,labels,num2str(j),'A1'); 
        xlswrite(Etdatafile,to_write,num2str(j),'A2'); 
        %Add these curves to the summary plot for this defect
        figure(allEk); 
        subplot(2,2,1); 
        semilogy(Et_now{j,1},k_now{j,1},'LineWidth',3,'Color',cm(j,:));
        hold all; 
        subplot(2,2,2);
        semilogy(Et_now{j,2},k_now{j,2},'LineWidth',3,'Color',cm(j,:));
        hold all; 
        subplot(2,2,3); 
        semilogy(Et_now{j,1},1./alphanN_now{j,1},'LineWidth',3,'Color',cm(j,:)); 
        hold all; 
        subplot(2,2,4);
        semilogy(Et_now{j,2},1./alphanN_now{j,2},'LineWidth',3,'Color',cm(j,:)); 
        hold all; 
        %Gather the values at the intrinsic energy for this defect
        %determine which defect is dominant
        simple_summary(j,:) = [interp1(Et_now{j,1},k_now{j,1},0),...
            interp1(Et_now{j,1},alphanN_now{j,1},0),...
            interp1(Et_now{j,2},k_now{j,2},0),...
            interp1(Et_now{j,2},alphanN_now{j,2},0)];
        if type == 'p'
            Xs = (n0+deltans)./(p0+deltans);
        elseif type == 'n'
            Xs = (p0+deltans)./(n0+deltans);
        end
        indexFZ = find(~cellfun(@isempty,defect_now(j,:)));
        if isempty(indexFZ)==0
            %assume we only have one FZ index, this could be a problem later
            all_defects = defect_now{j,indexFZ}; 
            full_details = all_defects{2}; 
            X = full_details{6}; tauSRH = full_details{7}; 
            for m = 1:length(Xs) 
                tauSRH_now = interp1(X,tauSRH,Xs(m)); 
                def1 = Xs(m)*fits_now(1,1)+fits_now(1,2); 
                def2 = Xs(m)*fits_now(1,3)+fits_now(1,4); 
                if abs(def1-tauSRH_now)<abs(def2-tauSRH_now)
                    dominance(j,m) = 1; 
                elseif abs(def2-tauSRH_now)<abs(def1-tauSRH_now)
                    dominance(j,m) = 2; 
                else
                    dominance(j,m) = 5; 
                end
            end
        else
            dominance(j,:) = 5; 
        end
    end
    
    %write the dominance and simple summary data for all the times
    xlswrite(Etdatafile,[simple_summary,dominance],num2str(0),'A2'); 
    labels = {'k defect 1','alphanN defect 1','k defect 2','alphanN defect 2'}; 
    for m = 1:length(deltans)
        labels{end+1} = ['dominant at ' num2str(deltans(m),'%0.e')]; 
    end
    xlswrite(Etdatafile,labels,num2str(0),'A1');   
    
    %Make a summary figure of k values etc with dominance indicated
    deltans_fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for m = 1:length(deltans)
        figure(deltans_fig); clf; 
        subplot(2,2,1); 
        semilogx(data(:,1),simple_summary(:,1),'rx-','LineWidth',2,'MarkerSize',8); 
        hold all; 
        subplot(2,2,2); 
        semilogx(data(:,1),1.\simple_summary(:,2),'rx-','LineWidth',2,'MarkerSize',8); 
        hold all; 
        subplot(2,2,3); 
        semilogx(data(:,1),simple_summary(:,3),'bx-','LineWidth',2,'MarkerSize',8); 
        hold all;
        subplot(2,2,4); 
        semilogx(data(:,1),1.\simple_summary(:,4),'bx-','LineWidth',2,'MarkerSize',8); 
        hold all;
        for h = 1:times
            if dominance(h,m) ~= 5 && dominance(h,m) == 1
                subplot(2,2,1); 
                semilogx(data(h,1),simple_summary(h,1),'ro','LineWidth',2,'MarkerSize',8); 
                hold all;
                subplot(2,2,2); 
                semilogx(data(h,1),1\simple_summary(h,2),'ro','LineWidth',2,'MarkerSize',8); 
                hold all; 
            elseif dominance(h,m) ~= 5 && dominance(h,m) == 2
                subplot(2,2,3); 
                semilogx(data(h,1),simple_summary(h,3),'bo','LineWidth',2,'MarkerSize',8); 
                hold all;
                subplot(2,2,4); 
                semilogx(data(h,1),1\simple_summary(h,4),'bo','LineWidth',2,'MarkerSize',8);
                hold all; 
            end
        end
        subplot(2,2,1); 
        xlabel('time [s]'); 
        ylabel('midgap k [-]'); 
        title('defect 1'); 
        subplot(2,2,2);
        xlabel('time [s]'); 
        ylabel('midgap \tau_n_0 [\mus]'); 
        title('defect 1');
        subplot(2,2,3); 
        xlabel('time [s]'); 
        ylabel('midgap k [-]'); 
        title('defect 2'); 
        subplot(2,2,4);
        xlabel('time [s]'); 
        ylabel('midgap \tau_n_0 [\mus]'); 
        title('defect 2');
        save_this = [samp_dir{2,i} '\Ektau_midgap_' num2str(deltans(m),'%0.e')];
        hgsave(deltans_fig,save_this);
        print(deltans_fig,'-dpng','-r0',[save_this '.png']);
    end
    
    figure(allEk); 
    subplot(2,2,1);
    xlabel('E_t-E_i [eV]'); 
    ylabel('k [-]'); 
    xlim([min(Et_now{1,1}) max(Et_now{1,1})]); 
    legend(num2str(data(:,1)),'Location','NorthEastOutside'); 
    title('default defect 1'); 
    subplot(2,2,2); 
    xlabel('E_t-E_i [eV]'); 
    ylabel('k [-]'); 
    xlim([min(Et_now{1,1}) max(Et_now{1,1})]);
    title('default defect 2');
    subplot(2,2,3); 
    xlabel('E_t-E_i [eV]'); 
    ylabel('\tau_{n0} [s]'); 
    xlim([min(Et_now{1,1}) max(Et_now{1,1})]);
    subplot(2,2,4); 
    xlabel('E_t-E_i [eV]'); 
    ylabel('\tau_{n0} [s]'); 
    xlim([min(Et_now{1,1}) max(Et_now{1,1})]);
    save_this = [samp_dir{2,i} '\fullEktau_over_time'];
    hgsave(gcf,save_this);
    print(gcf,'-dpng','-r0',[save_this '.png']);
    
    %Make a summary of fit parameters over time for each defect, including
    %a histogram
    fit_time=figure('units','normalized','outerposition',[0 0 1 1]); 
    subplot(4,2,1); 
    semilogx(data(:,1),data(:,2),'ro-','LineWidth',2); 
    xlabel('time [s]'); 
    ylabel('m'); 
    title('defect 1 slope'); 
    subplot(4,2,3); 
    semilogx(data(:,1),data(:,3),'ro-','LineWidth',2); 
    xlabel('time [s]'); 
    ylabel('b'); 
    title('defect 1 intercept'); 
    subplot(4,2,5); 
    semilogx(data(:,1),data(:,4),'bo-','LineWidth',2); 
    xlabel('time [s]'); 
    ylabel('m'); 
    title('defect 2 slope'); 
    subplot(4,2,7); 
    semilogx(data(:,1),data(:,5),'bo-','LineWidth',2); 
    xlabel('time [s]'); 
    ylabel('b'); 
    title('defect 2 intercept'); 
    subplot(4,2,2); 
    hist(data(:,2),nbins); 
    h = findobj(gca,'Type','patch');
    h.FaceColor = 'r'; h.EdgeColor = 'w';
    subplot(4,2,4); 
    hist(data(:,3),nbins); 
    h = findobj(gca,'Type','patch');
    h.FaceColor = 'r'; h.EdgeColor = 'w';
    subplot(4,2,6); 
    hist(data(:,4),nbins); 
    h = findobj(gca,'Type','patch');
    h.FaceColor = 'b'; h.EdgeColor = 'w';
    subplot(4,2,8); 
    hist(data(:,5),nbins); 
    h = findobj(gca,'Type','patch');
    h.FaceColor = 'b'; h.EdgeColor = 'w';
    save_this = [samp_dir{2,i} '\fits_over_time'];
    hgsave(fit_time,save_this);
    print(fit_time,'-dpng','-r0',[save_this '.png']);
    
    close all; 
end