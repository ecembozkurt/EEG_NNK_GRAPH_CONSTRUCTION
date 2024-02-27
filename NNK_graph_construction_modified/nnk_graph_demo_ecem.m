%clc;
%clear;
%close all;
%% read data
%knn_param=496; k_choice=496; reg=1e-6;
%X=newA;
%X=eucdis;

%X=degsim;
%X=coeff;
N = size(X,2);


%% compute graphs
tic
D = DistEuclideanPiotrDollar(X,X); % pairwise squared Euclidean distances

D(find(abs(D)<reg))=0; %ecem
D=abs(D);

%D=cossim-diag(ones(N,1));

directed_knn_mask = sparse(GD_BuildDirectedKnnGraph(D,knn_param,'dist'));
distance_mask_time = toc;
[kD,ix_kD] = sort(abs(sqrt(D)), 'ascend');
sigma2 = full(mean(kD(k_choice, :)));% for in 3 sigma
%sigma2=1; %doesn't dep on K
%*1e-1
G = exp(-D./(2*sigma2*sigma2));


%G=D;
% aecem=1;
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

similarity_time = toc;
knn_mask = max(directed_knn_mask, directed_knn_mask');
%knn_mask=directed_knn_mask; %ecem's trial, directed
symmetrization_time = toc;
%%
fprintf('Computing the adj and L of %d-NN graph with sigma %0.4f...\n', knn_param, sigma2)
W_knn = G .* knn_mask;
%W_knn=knn_mask;
W_knn(W_knn<reg) = 0;
knn_time = toc;

%%
fprintf('Computing the adj and L NNK...\n')
tic
condi=cond(G);
G2=log(1./D); 
G2(find(G2==Inf))=0;
G2=G2./max(max(G2));
G2=exp(-D);
if ~isnan(G)
    [W_nnk,error] = nnk_inverse_kernel_graph(G, directed_knn_mask, knn_param, reg); % choose the min k-NN sim

    %[W_nnk,error] = nnk_inverse_kernel_graph(G,ones(size(G,1),size(G,2)), size(G,1)-1, reg); % choose the min k-NN sim

else
    W_nnk=0;
end
nnk_time = toc + similarity_time;

%%
time_values = {knn_time, nnk_time};
sparsity_values = {length(find(W_knn))/2, length(find(W_nnk))/2};

W_nnk_cossim=W_nnk;
W_knn_cossim=W_knn;

% [rel1,rel2,edge]=find(W_nnk);
% G=graph(rel1,rel2,edge);
% LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
% figure(99),plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)
