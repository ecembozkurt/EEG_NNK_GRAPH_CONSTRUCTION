%clc;
%clear;
%close all;
%% read data
%knn_param=496; k_choice=496; reg=1e-6;
X=cossim;
sum(find(X>0.99));
N = size(X,2);
%results_folder = ['results3/']; 
%dir_result = mkdir(results_folder);

D = DistEuclideanPiotrDollar(X,X); % pairwise squared Euclidean distances
%condg=cond(D)
D(find(D<0))=0; %ecem
directed_knn_mask = sparse(GD_BuildDirectedKnnGraph(D,knn_param,'dist'));
D_sort=sort(D, 'ascend');
kD = sort(sqrt(D), 'ascend');
%figure(66666),plot(kD),hold all,
sigma = full(mean(kD(k_choice, :)));%/3;
G = exp(-D./(2*sigma*sigma));
%figure,imagesc(G),title('Matrix G'),colorbar
%local scaling kernel
aecem=1;
condi=Inf;
while condi>10
    
for i=1:size(X,1)
    for j= 1:size(X,2)
        sigma_i=kD(7,i);
        sigma_j=kD(7,j);
        G(i,j)=exp(-D(i,j)./(aecem*sigma_i*sigma_j));
        
    end
end

condi=cond(G)
aecem=aecem/100;
end

