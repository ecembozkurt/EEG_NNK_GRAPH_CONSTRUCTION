%node_stats:
W_NNK_ALL_R=W_M;
W_NNK_ALL_M=W_M;

%degree per node

%distribution of sensors for each

figure(888),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    imagesc(squeeze(W_NNK_ALL_R(:,i,:))'),hold on,
    title(sprintf('Channel %d',i))
end
sgtitle('REST: Weights for each channel')
hold off;

figure(777),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    imagesc(squeeze(W_NNK_ALL_M(:,i,:))'),hold on,
    title(sprintf('Channel %d',i))
end
sgtitle('MUSIC: Weights for each channel')
hold off;

SUM_weights_R=sum(W_NNK_ALL_R,1);
SUM_weights_M=sum(W_NNK_ALL_M,1);

mean_R=(SUM_weights_R./T);
mean_M=(SUM_weights_M./T);


var_R=zeros(ch_number,ch_number);
var_M=zeros(ch_number,ch_number);
sec_min_R=zeros(ch_number,ch_number);
sec_min_M=zeros(ch_number,ch_number);
median_R=zeros(ch_number,ch_number);
median_M=zeros(ch_number,ch_number);
sec_max_R=zeros(ch_number,ch_number);
sec_max_M=zeros(ch_number,ch_number);
mean_conn_R=zeros(ch_number,ch_number);
mean_conn_M=zeros(ch_number,ch_number);
for i=1:ch_number
    for j=1:ch_number
        var_R(i,j)=var(W_NNK_ALL_R(:,i,j));
        var_M(i,j)=var(W_NNK_ALL_M(:,i,j));
        
        if i~=j
            sort_R = sort(nonzeros(W_NNK_ALL_R(:,i,j)));
            sort_M = sort(nonzeros(W_NNK_ALL_M(:,i,j)));
            if ~isempty(sort_R)
                sec_min_R(i,j) = sort_R(1);
                sec_max_R(i,j) = sort_R(end);
                
                sec_min_M(i,j) = sort_M(1);
                sec_max_M(i,j) = sort_M(end);
            end
        end
        mean_conn_R(i,j)=mean(W_NNK_ALL_R(:,i,j)>0);
        mean_conn_M(i,j)=mean(W_NNK_ALL_M(:,i,j)>0);
        
        median_R(i,j)=median(W_NNK_ALL_R(:,i,j));
        median_M(i,j)=median(W_NNK_ALL_M(:,i,j));
    end
end

figure(444),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    plot(var_R(i,:)),hold on,
    plot(var_M(i,:))
    title(sprintf('Channel %d',i))
    legend('Rest','Music')
end
sgtitle('Variance')



figure(111),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    plot(squeeze(mean_R(1,i,:))),hold on,
    plot(squeeze (mean_M(1,i,:)))
    title(sprintf('Channel %d',i))
    legend('Rest','Music')
end
sgtitle('Mean')

figure(333),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    plot((sec_min_R(i,:))),hold on,
    plot((sec_min_M(i,:)))
    title(sprintf('Channel %d',i))
    legend('Rest','Music')
end
sgtitle('Nonzero Min')


figure(222),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    plot(median_R(i,:)),hold on,
    plot(median_M(i,:))
    title(sprintf('Channel %d',i))
    legend('Rest','Music')
end
sgtitle('Median')

figure(333),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    plot(sec_max_R(i,:)),hold on,
    plot( sec_max_M(i,:))
    title(sprintf('Channel %d',i))
    legend('Rest','Music')
end
sgtitle('Second Max')

figure(999),
for i=1:ch_number
    subplot(ceil(ch_number/2),2,i)
    plot(mean_conn_R(i,:)),hold on,
    plot(mean_conn_M(i,:))
    title(sprintf('Channel %d',i))
    legend('Rest','Music')
end
sgtitle('Mean Connectivity')

%MSE
figure(666),
for i=1:ch_number
    sb(i)=subplot(ceil(ch_number/2),2,i);
    plot(sqrt(squeeze(mean_R(1,i,:)-mean_M(1,i,:)).^2))
    title(sprintf('Channel %d',i))
    %legend('Rest','Music')
    grid on
end
sgtitle('MSE of Mean')
linkaxes(sb,'y')
ax = axis;
axis([ax(1:2)  0  0.018])

