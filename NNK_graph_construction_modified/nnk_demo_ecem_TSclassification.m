X=SimMatrix;
X=(X-min(min(X)))./(max(max(X))-min(min(X)));
X(isnan(X))=1;

N = size(X,2);
reg=1e-6;

%% compute graphs
D = DistEuclideanPiotrDollar(X,X); % pairwise squared Euclidean distances
cd '../Kernel-Sparse-Representation-master/KernelKSVD for time series';%%DWT, ERP,TWED
%D = dist_DTW(X,X,1);
%D= dist_ERP(X,X,1);
%D= dist_TWED(X,X,0.1,0);
cd '../../NNK_graph_construction-master'
D(find(abs(D)<reg))=0; %ecem
D=abs(D);
knn_param=1;
%D=cossim-diag(ones(N,1));
figure,imagesc(D),title('Gaussian Kernel: DTW')
CCC=0;
old_W_nnk=zeros(N,N);
W_nnk=zeros(N,N);
knn_par= ch-1;
while CCC<1e3 
    old_W_nnk=W_nnk;
    directed_knn_mask = sparse(GD_BuildDirectedKnnGraph(D,knn_par,'dist'));
    [kD,ix_kD] = sort(abs(sqrt(D)), 'ascend');
    sigma2 = full(mean(kD(knn_par+1, :)));% for in 3 sigma
    %sigma2=1; %doesn't dep on K
    %*1e-1
    G = exp(-D./(2*sigma2*sigma2));
    CCC=cond(G);
    
    %G=D;
    % aecem=1;35
    % condi=Inf;
    % while condi>1e4
    %
    % for i=1:size(X,1)
    %     for j= 1:size(X,2)
    %         sigma_i=kD(k_choice,i);
    %         sigma_j=kD(k_choice,j);
    %         G(i,j)=exp(-D(i,j)./(aecem*sigma_i*sigma_j));
    %
    %     end
    % end
    %
    % condi=cond(G)
    % aecem=aecem/100;
    % end
    
    knn_mask = max(directed_knn_mask, directed_knn_mask');
    %knn_mask=directed_knn_mask; %ecem's trial, directed
    %%
%    fprintf('Computing the adj and L of %d-NN graph with sigma %0.4f...\n', knn_param, sigma2)
    W_knn = G .* knn_mask;
    W_knn(W_knn<reg) = 0;
    
    %%
 %   fprintf('Computing the adj and L NNK...\n')
    
    %G2=log(1./D);
    %G2(find(G2==Inf))=0;
    %G2=G2./max(max(G2));
   % G2=exp(-D);
    if ~isnan(G)
        [W_nnk,error] = nnk_inverse_kernel_graph(G, directed_knn_mask, knn_par, reg); % choose the min k-NN sim
        %[W_nnk,error,ERR_ecem] = omp_nnk_inverse_kernel_graph(G,directed_knn_mask,knn_par,reg);
        %[W_nnk,error] = omp_nnk_inverse_kernel_graph_group(G,directed_knn_mask,knn_par,reg,true,maxdelay); 
       
    else
        W_nnk=0;
    end

    knn_param=knn_param+1;
    disp(knn_param)

    if W_nnk==0
        W_nnk=old_W_nnk;
    end
    if length(find(old_W_nnk))== length(find(W_nnk))
        length(find(W_nnk))
        break
    end
    W_nnk=W_nnk;
end

% figure,
% for iii=1:ch_number
%     subplot(7,2,iii),plot(ERR_ecem(iii,:),'LineWidth',1),xlabel('k'),ylabel('Error')
%     title(sprintf('Node: %d',iii))
% end
%%
disp('knn param')
    disp(knn_param)
sparsity_values = {length(find(W_knn)), length(find(W_nnk))};

W_nnk_cossim=max(W_nnk,W_nnk');
W_knn_cossim=max(W_knn,W_knn');
