% MAIN DEMO for EEG Graph Construction with NNK Alogrithm
% Author: Ecem
% Date: 02/2024 (updated)
% Input: EEG Signals
% Output: Graph weights (W_NNK) 
%=================================================================

%% LOAD DATA:
node_locations=load('../data/eeg_xyz.mat');
load('../data/music04_alpha.mat')
music_data=data;
load('../data/rest04_alpha.mat')
rest_data=data;
ch_number=14;
time_window=[500];%100,500,1000,
slide_step=[100];
addpath('../../NNK_graph_construction-master/')
for RM=1:2
    if RM==1
        A=rest_data;
        sData='REST ';
    end
    
    if RM==2
        A=music_data;
        sData='MUSIC ';
    end
    [ch,T]=size(A);
    
    
    %% INITIALIZE VARIABLES:
    frame_number=floor((T-time_window(1))/(slide_step));
    SIM=zeros(ch,ch, frame_number);
    prev_W_nnk=zeros(ch_number,ch_number);
    W_nnk=zeros(ch_number,ch_number);
    eeg_xyz=load('../data/emotiv14_xyz.mat').emotiv14_xyz;

    %% Geometric Relation
    distances=zeros(ch_number,ch_number);
    for i=1:ch_number
        for j=1:ch_number
            distances(i,j)=euclidianDist_3D(eeg_xyz(i,:),eeg_xyz(j,:));
        end
    end
    
    
    GaussKernel= GaussianKernel(distances,500);

    for i =1:length(time_window)  % for each time window size
        %% Number of Connections
        
        N_CONN=zeros(frame_number,ch_number);
        SUMW=zeros(frame_number,ch_number);

         %% SEGMENTATION:
        for k=1:frame_number % slide window
            newA=A(1:ch_number,(1+(k-1)*slide_step):(time_window(i)+slide_step*(k-1)));
            
            %% FIND SIMILARITY : COS SIM
            cossim=channelCosSim(newA);
            %SimMatrix=getMI(newA);
            %SimMatrix=getPLV(newA)
 
            SimMatrix=cossim;
            SIM(:,:,k)=SimMatrix;
            
            
            %% NNK GRAPH CONSTRUCTION:
            prev_W_nnk=W_nnk;
            
           [W_nnk, W_knn]=nnk_graph_construction_ecem(SimMatrix, 10, 10, 1e-6);

           % cd ('../SPIS-Resting-State-Dataset-master-2/')
       
            [sumW,n_conn]=plotChange_timeseries(W_nnk);
            N_CONN(k,:)=n_conn;
            SUMW(k,:)=sumW;
        end
        
        
        %% VISUALIZE: PLOTS

        fig=figure;
        imagesc(N_CONN')
        colorbar;
        stitle=sprintf('Time Window: %d, Step Size %d',time_window(i),slide_step);
        title(sData,stitle)
        xlabel('index')
        ylabel('Channels')
        set(fig,'Position',[1000 0 1200 400])
        fname=sprintf('./EEGCOS%s_%s.jpg',sData,stitle);
        saveas(fig,fname)

        
        %fig2=figure;
        %imagesc(SimMatrix)
        %colorbar;
        %stitle=sprintf('EUC: Time Window: %d, Step Size %d',time_window(i),slide_step);
        %sEuc=sprintf('Euclidian Distances');
        %title(sData,stitle)
        %xlabel('Channels')
        %ylabel('Channels')
        %fname2=sprintf('./EEGCOS%s_%s.jpg',sData,stitle);
        %saveas(fig2,fname2)
    en                                                      d
    
    figure(700+RM),imagesc(reshape(SIM,ch*ch,k))
    title('SIM Matrix through time frames')
    xlabel('frame number')
    ylabel('ch-ch rel.')
    
end


%%  METRICS AND STATISTICAL ANALYSIS OF RESULTS:














