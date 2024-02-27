cd ('../EEG-W-NNK-main/')
%[V,a1,c1,a2,c2,a3,c3,ph1,ph2,ph3]=gen_syntheticdata(20,1000);
[V,a1,c1,a2,c2,a3,c3,ph1,ph2,ph3]=gen_syntheticdata_phase(20,1000);

time_window=100;
slide_step=time_window;
condi=0;
T1=size(V,2);
ch_number=size(V,1);
frame_number=floor((T1-time_window)/(slide_step));
if frame_number==0
    frame_number=1;
end
frame_number=1;
SIM=zeros(ch_number,ch_number, frame_number);
prev_W_nnk=zeros(ch_number,ch_number);
W_nnk=zeros(ch_number,ch_number);

%% Number of Connections
N_CONN=zeros(frame_number,ch_number);
SUMW=zeros(frame_number,ch_number);
DECONS=zeros(frame_number,1);
DECONS2=zeros(frame_number,1);
W_M=zeros(frame_number,ch_number,ch_number);
axs=zeros(frame_number,time_window);
W_M_stitch=zeros(T1,ch_number,ch_number);
stitch_M=zeros(frame_number,T1);
newA=zeros(ch_number,time_window);
PREVA=zeros(frame_number,ch_number,time_window);

tw_inx=0; knn_inx=0;
for time_window=time_window%[10,20,50,100,200,500]
    tw_inx=tw_inx+1;
    for knnp=10%[10,20,50,100,300,400,450,469]
        knn_inx=knn_inx+1;
        for k=1:3% slide window
            prevA=newA; %prev. wind
            fr_indices=(1+(k-1)*slide_step):(time_window+slide_step*(k-1));
            newA=V(1:ch_number,fr_indices);
            for  iix=1: time_window
                VG_2{iix}=fast_NVG(newA(:,iix),1:ch_number,'w',0);
            end
            for ii=1:ch_number
                
                VG{ii}= fast_NVG(newA(ii,:),fr_indices,'w',0);
                %VG{ii}= fast_HVG(newA(ii,:),fr_indices,0);
            end
            
            PREVA(k,:,:)=newA;
            if sum(sum(isnan(newA)))>0
                disp('NaN error')
            end
            
            %% COS SIM
            cd ('../EEG-W-NNK-main/')
       
            degsim=channelDeg2d(VG);
            %cossim=channelCosSim2d(VG);
            cossim=channelCosSim(degsim);
            for i4=1:ch_number
                for i5=1:ch_number
                    [~,lag(i4,i5)]=synch_DS(degsim(i4,:),degsim(i5,:),10);
                end
            end
            cossimL=channelCosSimDelay(degsim,lag);
            SimMatrix=degsim;
            SIM=degsim;
            cd ('../NNK_graph_construction-master/')
            %% NNK
            prev_W_nnk=W_nnk;
            %%Ecem'S
            knn_param=knnp; k_choice=knnp; reg=1e-6;
            nnk_demo_ecem_incremental
            
            cd ('../EEG-W-NNK-main/')
            %%PLOT
            [rel1,rel2,edge]=find(W_nnk);
            Graph0=graph(rel1,rel2,edge);
            LWidths = 5*Graph0.Edges.Weight/max(Graph0.Edges.Weight);
            figure(200+k),hh=plot(Graph0,'LineWidth',LWidths) %'EdgeLabel',G.Edges.Weight,
            title('Proposed NNK Method')
            
            for idx=1:max(max(rel1),max(rel2))
                switch c1(idx)/0.05
                    case 0.05
                        highlight(hh,idx,'NodeColor','green')
                    case 0.1
                        highlight(hh,idx,'NodeColor','black')
                    case 0.15
                        highlight(hh,idx,'NodeColor','red')
                    case 0.2
                        highlight(hh,idx,'NodeColor','yellow')
                    case 0.25
                        highlight(hh,idx,'NodeColor',[0.1 0.5 0.5])
                    case 0.3
                        highlight(hh,idx,'NodeColor','cyan')
                    case 0.35
                        highlight(hh,idx,'NodeColor',[0 0.5 1])
                    case 0.4
                        highlight(hh,idx,'NodeColor',[1 0.5 0.5])
                    case 0.45
                        highlight(hh,idx,'NodeColor',[1 0.5 1])
                    case 0.5
                        highlight(hh,idx,'NodeColor',[0.5 0 1])
                    case 0.55
                        highlight(hh,idx,'NodeColor',[0.5 1 1])
                    case 0.6
                        highlight(hh,idx,'NodeColor',[0.5 1 0.5])
                    case 0.65
                        highlight(hh,idx,'NodeColor',[0.25 1 1])
                    case 0.7
                        highlight(hh,idx,'NodeColor',[0 0.25 1])
                    case 0.75
                        highlight(hh,idx,'NodeColor',[1 0.25 0.25])
                    case 0.8
                        highlight(hh,idx,'NodeColor',[1 0.25 1])
                    case 0.85
                        highlight(hh,idx,'NodeColor',[0.25 0 1])
                    case 0.9
                        highlight(hh,idx,'NodeColor',[0.7 0.7 0.7])
                    case 01
                        highlight(hh,idx,'NodeColor',[0.75 0.75 1])
                    otherwise
                        disp(idx)
                end
            end
            
            
            %%
            
            for iiii=7:15%ch_number
                
                [rel3,rel4,edge3]=find(abs(full(VG{iiii})));
                Graph02=graph(rel3,rel4,edge3);
                LWidths = 5*Graph02.Edges.Weight/max(Graph02.Edges.Weight);
                %figure(9900+k),subplot(10,2,2*iiii-13),h2=plot(Graph02,'LineWidth',LWidths) %'EdgeLabel',G.Edges.Weight,
                figure(9900+k),subplot(10,2,2*iiii-13),imagesc(abs(full(VG{iiii}))) %'EdgeLabel',G.Edges.Weight,
                title(sprintf("Node: %d ",iiii))
                figure(9900+k),subplot(10,2,2*iiii-12), plot(newA(iiii,:))
                title(sprintf("Node: %d ",iiii))
                sgtitle('VG')
                
            end
            %%
            [rel5,rel6,edge3]=find(W_knn);
            Graph03=graph(rel5,rel6,edge3);
            LWidths = 5*Graph03.Edges.Weight/max(Graph03.Edges.Weight);
            figure(400+k),hh3=plot(Graph03,'LineWidth',LWidths) %'EdgeLabel',G.Edges.Weight,
            title('KNN Method')
            
            for idx=1:max(max(rel5),max(rel6))
                switch c1(idx)/0.05
                    case 0.05
                        highlight(hh3,idx,'NodeColor','green')
                    case 0.1
                        highlight(hh3,idx,'NodeColor','black')
                    case 0.15
                        highlight(hh3,idx,'NodeColor','red')
                    case 0.2
                        highlight(hh3,idx,'NodeColor','yellow')
                    case 0.25
                        highlight(hh3,idx,'NodeColor',[0.1 0.5 0.5])
                    case 0.3
                        highlight(hh3,idx,'NodeColor','cyan')
                    case 0.35
                        highlight(hh3,idx,'NodeColor',[0 0.5 1])
                    case 0.4
                        highlight(hh3,idx,'NodeColor',[1 0.5 0.5])
                    case 0.45
                        highlight(hh3,idx,'NodeColor',[1 0.5 1])
                    case 0.5
                        highlight(hh3,idx,'NodeColor',[0.5 0 1])
                    case 0.55
                        highlight(hh3,idx,'NodeColor',[0.5 1 1])
                    case 0.6
                        highlight(hh3,idx,'NodeColor',[0.5 1 0.5])
                    case 0.65
                        highlight(hh3,idx,'NodeColor',[0.25 1 1])
                    case 0.7
                        highlight(hh3,idx,'NodeColor',[0 0.25 1])
                    case 0.75
                        highlight(hh3,idx,'NodeColor',[1 0.25 0.25])
                    case 0.8
                        highlight(hh3,idx,'NodeColor',[1 0.25 1])
                    case 0.85
                        highlight(hh3,idx,'NodeColor',[0.25 0 1])
                    case 0.9
                        highlight(hh3,idx,'NodeColor',[0.7 0.7 0.7])
                    case 01
                        highlight(hh3,idx,'NodeColor',[0.75 0.75 1])
                    otherwise
                        disp(idx)
                end
            end
            
            %%
        end
    end
end
