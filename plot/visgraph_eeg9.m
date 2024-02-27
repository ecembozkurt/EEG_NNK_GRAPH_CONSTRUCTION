eeg_xyz=load('../data/eeg_xyz_1020.mat').eeg_xyz1020;
load('../data/EEG9data.mat')
load('../data/interp_annot.mat')
ch_number=19;
addpath '../../Neonatal_Seizure_Detection-master/Neonatal_Seizure_Detection_Algorithm/neonatal_sez_det'

%% Preprocess
prep_eeg=zeros(size(EEG9));
for i=1:ch_number
    prep_eeg(i,:)=preprocess(EEG9(i,:));
end
A=prep_eeg;
A=prep_eeg(:,1:256:size(A,2));
T=length(A(1,:));

time_window=[50];%100,500,1000,
slide_step=[50];
frame_number=floor((T-time_window(1))/(slide_step))+1;
%SIM=zeros(ch_number,ch_number, frame_number);
prev_W_nnk=zeros(ch_number,ch_number);
W_nnk=zeros(ch_number,ch_number);
W_NNK=zeros(frame_number,ch_number,ch_number);
%% Geometric Relation
distances=zeros(ch_number,ch_number);
for i=1:ch_number
    for j=1:ch_number
        distances(i,j)=euclidianDist_3D(eeg_xyz(i,:),eeg_xyz(j,:));
    end
end


%GaussKernel= GaussianKernel(distances,500);

%figure, imagesc(distances)
%title('Euclidian Distances (mm)')
%figure, imagesc(GaussKernel)
%title('Gaussian Kernel')
for i =1:length(time_window)  % for each time window size
    %% Number of Connections
    
    N_CONN=zeros(frame_number,ch_number);
    % SUMW=zeros(frame_number,ch_number);
    %DECONS=zeros(frame_number,1);
    DECONS2=zeros(frame_number,1);
    
    for k=1:frame_number % slide window
        fr_indices=(1+(k-1)*slide_step):(time_window(i)+slide_step*(k-1));
        for ii=1:ch_number
            EEG_VG{ii}= fast_NVG(A(ii,fr_indices),fr_indices,'u',0);
        end
        
        %k;
        % (k-1)*time_window(i)
        % k*time_window(i)
        %[(1+(k-1)*slide_step),(time_window(i)+slide_step*(k-1))];
        %newA=A(1:ch_number,(1+(k-1)*slide_step):(time_window(i)+slide_step*(k-1)));
        annotate_wind=max(inrped_annot(1,(1+(k-1)*slide_step):(time_window(i)+slide_step*(k-1))));
        ANN(k)=annotate_wind;
        %% COS SIM
        cossim= channelCosSim2d(EEG_VG);
        %eucdis= channelEucDist2d(EEG_VG);
        %cossim=channelCosSim(newA);
        %EucDis0=getEuclidian(newA)./500;
        %EucDis=(ones(ch_number,ch_number)./(ones(ch_number,ch_number)+EucDis0));
        %        BD=getBD(newA);
        %SimMatrix=getMI(newA);
        %SimMatrix=getPLV(newA);
        
        %SimMatrix=EucDis;
        %SimMatrix=abs(BD);
        %SimMatrix=cossim;
        % figure, imagesc(cossim)
        %SIM(:,:,k)=SimMatrix;
        cd ('../NNK_graph_construction-master/')
        %% NNK
        prev_W_nnk=W_nnk;
        % [Y,W_lle] = lle(newA,19, 18);
        knn_param=18; k_choice=18; reg=1e-6;
        nnk_graph_demo_ecem
        cd ('../EEG-W-NNK-main/')
        W_NNK(k,:,:)=full(W_nnk);
        %title=sprintf('Time Window %d, Time Slot  between %d and %d', time_window(i),((k-1)*time_window(i)+1),k*time_window(i));
        %fig=plotGraph64(W_nnk, title);
        %
        %W_nnk=W_lle;
        [sumW,n_conn]=plotChange_timeseries(W_nnk, 'number of connections');
        N_CONN(k,:)=n_conn;
        %SUMW(k,:)=sumW;
        %%   [decon,~,~]=deltacon0(W_nnk,prev_W_nnk);
        %%    DECONS(k)=decon;
        %% plot
%         [rel1,rel2,edge]=find(W_nnk);
%         Graph0=graph(rel1,rel2,edge);
%         LWidths = 5*Graph0.Edges.Weight/max(Graph0.Edges.Weight);
%         figure(200+k),hh=plot(Graph0,'LineWidth',LWidths) %'EdgeLabel',G.Edges.Weight,
%         title('Proposed NNK Method')
        %%
        
        cd ('../deltacon/')
        [decon2,simSTD, time, timeSTD]=DeltaCon('naive',W_nnk,prev_W_nnk,0.1);
        cd ('../EEG-W-NNK-main/')
        DECONS2(k)=decon2;
        %fig=plotGraphALL2(W_nnk,'',ch_number);
        %alltogether
        %         fig=plotGraphALL(W_nnk, 'ALL connections',ch_number);
        %         fname=sprintf('./Images2/newGall_fig_tw_%d_btw_%d_and_%d.jpg',time_window(i),((k-1)*time_window(i)+1),k*time_window(i));
        %         saveas(fig,fname)
        %         close
        % if (k>=2)
        % if (abs(n_conn-N_CONN(k-1,:))>0)
        % fig=plotGraph64(W_nnk,' ',ch_number);
        % fname=sprintf('./ImageM/Music/fig_tw_%d_btw_%d_and_%d.jpg',time_window(i),((k-1)*time_window(i)+1),k*time_window(i));
        % saveas(fig,fname)
        % close
        % end
        % end
        
    end
end
% use EEG_VG for NNK graphs
interp_decon=interp(DECONS2,fix(3550/size(DECONS2,1)));
figure,subplot(211),plot(interp_decon>0.56)
subplot(212),plot(max(annotat)')
sgtitle('DECON(thr=0.56)')
%
maxPred=squeeze(max(max(W_NNK,[],3),[],2));
interp_max_pred=interp(maxPred,fix(3550/size(maxPred,1)));
figure,subplot(211),plot(interp_max_pred>0.66)
subplot(212),plot(max(annotat)')
sgtitle('Max of Weights (thr =0.66)')

sum_pred=squeeze(sum(sum(W_NNK,3),2));
interp_sum_pred=interp(sum_pred,fix(3550/size(maxPred,1)));
figure,subplot(211),plot(interp_sum_pred>4.5)
subplot(212),plot(max(annotat)')
sgtitle('Sum of Weights (thr=4.5)')

interp_nconn=interp(sum(N_CONN,2),fix(3550/size(maxPred,1)));
figure,subplot(211),plot(interp_nconn>80)
subplot(212),plot(max(annotat)')
sgtitle('Number of Connections (thr= 80)')


figure,subplot(211),plot(interp_max_pred>0.66),hold all,
plot(interp_sum_pred>4.5),plot(interp_nconn>80),plot(interp_max_pred>0.66),plot(interp_decon>0.56)
subplot(212),plot(max(annotat)')

figure,subplot(10,2,20),plot(max(annotat)')
for i=1:ch_number
interp_W=interp(squeeze(sum(W_NNK(:,i,:),3)),fix(3550/size(maxPred,1)));
subplot(10,2,i),plot(interp_W)
end
sgtitle('Sum of Weights per Channel')

%clustering coefficient

%network entropy

%modularity

%degree distribution index
