
cd /Users/MYQUEEN/Downloads/EEG-W-NNK-main

[A,a1,c1,a2,c2,a3,c3,ph1,ph2,ph3]=gen_syntheticdata_phase(20,1000);

time_window=100;
slide_step=time_window;
condi=0;
T1=size(A,2);

ch_number=size(A,1);

%A=A+0.005*ones(ch_number,366);
T=length(A(1,:));
%time_window=[10];
%slide_step=10;%[10];
frame_number=floor((T-time_window(1))/(slide_step))+1;
maxdelay=2;
fs=1;
prev_W_nnk=zeros((maxdelay+1)*ch_number,(maxdelay+1)*ch_number);
W_nnk=zeros((maxdelay+1)*ch_number,(maxdelay+1)*ch_number);

%% Geometric Relation
distances=zeros(ch_number,ch_number);


for i =1%length(time_window)  % for each time window size
    %% Number of Connections
    N_CONN=zeros(frame_number,(maxdelay+1)*ch_number);
    DECONS2=zeros(frame_number,1);
    
    for k=5%frame_number%:frame_number % slide window
        fr_indices=(1+(k-1)*slide_step):(time_window(i)+slide_step*(k-1));
        
        for ii=1:ch_number
            newA(ii,:)=A(ii,fr_indices);
            SYNT_VG{ii}= fast_NVG(A(ii,fr_indices),fr_indices,'w',0);
            
        end
        %% COS SIM
        cd ('../EEG-W-NNK-main/')
        SYNT_VG_new={};sVGtitle={};shift=[];
        for iml= 0:(maxdelay)
            for ich=1:ch_number
                %TEMP_VG_new{end+1}= circshift(TEMP_VG{ich},[iml*fs,iml*fs]);
                %fr_indices+iml
                SYNT_VG_new{end+1}=fast_NVG(A(ich,fr_indices+iml),fr_indices+iml,'w',0);
                sVGtitle{end+1}=sprintf('node %d,  shift %d',ich,iml);
                shift(end+1)=iml;
            end
        end
        degsim=channelDeg2d(SYNT_VG_new);
        coeff=channelCLUSCOEFF2d(SYNT_VG_new);
        %cossim=channelCosSim2d(VG);
        %cossim=channelCosSim(coeff);
        SimMatrix=degsim;
        SIM=degsim;
        
        %% NNK
        prev_W_nnk=W_nnk;
        
        %% COS SIM
        
        % cossim= channelCosSim2d);
        cd ('../NNK_graph_construction-master/')
        knn_param=3; k_choice=3; reg=1e-6;
        nnk_demo_ecem_incremental
        
        W_NNK(k,:,:)=full(W_nnk);
        
        W_nnk_pruned=pruneDelay(W_nnk,maxdelay);
        W_knn_pruned=pruneDelay(W_knn,maxdelay);
        sparsity_values3 = {length(find(W_knn)), length(find(W_nnk))}
        sparsity_values4 = {length(find(W_knn_pruned)), length(find(W_nnk_pruned))}
        
        cd ('../EEG-W-NNK-main/')
        [sumW,n_conn]=plotChange_timeseries(W_nnk, 'number of connections');
        N_CONN(k,:)=n_conn;
        
        cd ('../deltacon/')
        [decon2,simSTD, time, timeSTD]=DeltaCon('naive',W_nnk,prev_W_nnk,0.1);
        cd ('../EEG-W-NNK-main/')
        DECONS2(k)=decon2;
        
        N_nnk(k)=interlayer_conn(W_nnk,ch_number);
        N_knn(k)=interlayer_conn(W_knn,ch_number);
        
        
        %%%
        cd '/Users/MYQUEEN/Downloads/EEG-W-NNK-main'
        SV3(k,:)=sparsity_values3; SV4(k,:)=sparsity_values4;
    end
    
    
end
%% NNK
fig=figure;
for iii=1:ch_number
    subplot(10,2,iii),plot(newA(iii,:),'LineWidth',2),title(sprintf('Node %d',iii))
end
fig=figure;
[X,Y,Z]=meshgrid(1:4,1:5,0:maxdelay);
p4=plot3(X(:),Y(:),Z(:),'*') ;

grid on
for i=1:size(W_nnk,1)
    hold all
    p4=plot3(rem(rem((i-1) ,ch_number),4)+1,floor(rem((i-1),ch_number)/4)+1,floor((i-1)/ch_number),'*');
    %[floor(rem(i-1,20)/5),rem(rem(i-1,20),4),floor(i-1/20)]
    text(rem(rem((i-1) ,ch_number),4)+1,floor(rem((i-1),ch_number)/4)+1,floor((i-1)/ch_number)+0.1,sVGtitle{i})
    
    for k=1:size(W_nnk,1)
        if W_nnk(i,k)==0 || W_nnk(i,k)==Inf
        else
            cc= abs(W_nnk(i,k));
            p4=plot3([rem(rem((i-1) ,ch_number),4)+1,rem(rem((k-1) ,ch_number),4)+1],[floor(rem((i-1),ch_number)/4)+1,floor(rem((k-1),ch_number)/4)+1],[floor((i-1)/ch_number),floor((k-1)/ch_number)]);
            grid on
            p4.LineWidth=cc*5;
        end
    end
end
title('NNK')
hold off;

%% KNN
fig=figure;
[X,Y,Z]=meshgrid(1:4,1:5,0:maxdelay);
p5=plot3(X(:),Y(:),Z(:),'*') ;
grid on
for i=1:size(W_knn,1)
    hold all
    p5=plot3(rem(rem((i-1) ,ch_number),4)+1,floor(rem((i-1),ch_number)/4)+1,floor((i-1)/ch_number),'*');
    %[floor(rem(i-1,20)/5),rem(rem(i-1,20),4),floor(i-1/20)]
    text(rem(rem((i-1) ,ch_number),4)+1,floor(rem((i-1),ch_number)/4)+1,floor((i-1)/ch_number)+0.1,sVGtitle{i})
    
    
    for k=1:size(W_knn,1)
        if W_knn(i,k)==0 || W_knn(i,k)==Inf
        else
            cc= abs(W_knn(i,k));
            p5=plot3([rem(rem((i-1) ,ch_number),4)+1,rem(rem((k-1) ,ch_number),4)+1],[floor(rem((i-1),ch_number)/4)+1,floor(rem((k-1),ch_number)/4)+1],[floor((i-1)/ch_number),floor((k-1)/ch_number)]);
            grid on
            p5.LineWidth=cc*5;
        end
    end
end
title('KNN')
hold off;

