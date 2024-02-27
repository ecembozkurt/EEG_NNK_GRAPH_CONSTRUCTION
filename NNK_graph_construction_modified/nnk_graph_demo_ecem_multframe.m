%clc;
%clear;
%close all;
%% read data
%knn_param=496; k_choice=496; reg=1e-6;
%X=newA;
X=squeeze(cossim(2,:,:));
Y=squeeze(cossim(1,:,:));
N = size(X,2);
%% compute graphs
D = DistEuclideanPiotrDollar(X,X); % pairwise squared Euclidean distances
D(find(abs(D)<reg))=0; %ecem
D=abs(D);

D2 = DistEuclideanPiotrDollar(Y,Y); % pairwise squared Euclidean distances
D2(find(abs(D2)<reg))=0; %ecem
D2=abs(D2);


directed_knn_mask = sparse(GD_BuildDirectedKnnGraph(D,knn_param,'dist'));
[kD,ix_kD] = sort(abs(sqrt(D)), 'ascend');
sigma2 = full(mean(kD(k_choice+1, :))); %/3 for in 3 sigma
%sigma2=1; %doesn't dep on K
G = exp(-D./(2*sigma2*sigma2*1e-2));
%%prev
directed_knn_mask2 = sparse(GD_BuildDirectedKnnGraph(D2,knn_param,'dist'));
[kD2,ix_kD2] = sort(abs(sqrt(D2)), 'ascend');
sigma2 = full(mean(kD2(10, :))); %/3 for in 3 sigma
%sigma2=1; %doesn't dep on K
G2 = exp(-D2./(2*sigma2*sigma2*1e-2));
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

knn_mask = max(directed_knn_mask, directed_knn_mask');
knn_mask2 = max(directed_knn_mask2, directed_knn_mask2');
%knn_mask=directed_knn_mask; %ecem's trial, directed
%%
fprintf('Computing the adj and L of %d-NN graph with sigma %0.4f...\n', knn_param, sigma2)
W_knn = G .* knn_mask;
W_knn2 = G2 .* knn_mask2;
%W_knn=knn_mask;
W_knn(W_knn<reg) = 0;
W_knn2(W_knn2<reg) = 0;

%%
fprintf('Computing the adj and L NNK...\n')

if ~isnan(G)
    [W_nnk,error] = nnk_inverse_kernel_graph(G, directed_knn_mask, knn_param, reg); % choose the min k-NN sim
else
    W_nnk=0; 
end


if ~isnan(G2)
     [W_nnk2,error2] = nnk_inverse_kernel_graph(G2, directed_knn_mask2, knn_param, reg); % choose the min k-NN sim
else
      W_nnk2=0;
end
    
%%


% [rel1,rel2,edge]=find(W_nnk);
% Grr=graph(rel1,rel2,edge);
% LWidths = 5*Grr.Edges.Weight/max(Grr.Edges.Weight);
% figure(99),plot(Grr,'EdgeLabel',Grr.Edges.Weight,'LineWidth',LWidths)
% 
% [rel3,rel4,edge2]=find(W_nnk2);
% Grr2=graph(rel3,rel4,edge2);
% LWidths = 5*Grr2.Edges.Weight/max(Grr2.Edges.Weight);
% figure(88),plot(Grr2,'EdgeLabel',Grr2.Edges.Weight,'LineWidth',LWidths)
