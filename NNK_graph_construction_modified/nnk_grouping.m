%grouping

X=cossim;
[i,j,val]=find(cossim>=0.999);

%
icc=0;
for in=1:length(i)
    %mask_group(i(in),j(in))=0;
    %W_nnk_group(i(in),j(in))=1;
    if(i(in)==j(in))
    else
       icc=icc+1;
       pick(icc,:)=[i(in),j(in)]; % similar weights pick
    end

end
pick(:,1)
pick(:,2)
X_pick_row=mean(X(pick(:,1),:),1);
X_pick_col=mean(X(:,pick(:,2)),2);
X(pick(:,1),:)=NaN;
X(:,pick(:,2))=NaN;
X_pick_row(pick(:,1))=[];
X_pick_col(pick(:,2))=[];
X_pick=(X_pick_row+X_pick_col')./2;
newX=X;
size(X)
newX(pick(:,1),:)=[];
newX(:,pick(:,2))=[];
size(newX)
% add the mean of the group to the end
newX(end+1,:)=X_pick;
newX(:,end+1)=[X_pick 1];
N=size(newX,1);


X=newX;




%%%%%%%%%%%%%%%%%%%%
%N = size(X,2);
mask_group=ones(N,N);

D = DistEuclideanPiotrDollar(X,X); % pairwise squared Euclidean distances
D(find(abs(D)<1e-11))=0; %ecem
%D=1-cossim;
if N<knn_param, knn_param=N; k_choice=N; end
directed_knn_mask = sparse(GD_BuildDirectedKnnGraph(D,knn_param,'dist'));
[kD,ix_kD] = sort(abs(sqrt(D)), 'ascend');
sigma2 = full(mean(kD(k_choice+1, :))); %/3 for in 3 sigma
%sigma2=1; %doesn't dep on K
G = exp(-D./(2*sigma2*sigma2));

D=abs(D);
threshold=1e-10;
%W_g=zeros(N,N);
W_g=1./(1+D);
newG=G;
%[i,j]=find(D<=threshold);

%directed_knn_mask_new=directed_knn_mask_new(1:N_gr,1:N_gr);

%figure,imagesc(mask_group)
%figure,imagesc(W_g)
%figure,imagesc(new_G)
directed_knn_mask_new= directed_knn_mask;
%W_nnk = nnk_inverse_kernel_graph(G, directed_knn_mask, knn_param, reg); % choose the min k-NN sim
W_nnk_group = nnk_inverse_kernel_graph(newG, directed_knn_mask_new, knn_param, reg); % choose the min k-NN sim
W_nnk_group=max(W_nnk_group,W_nnk_group');
W_nnk_group=W_nnk_group-eye(N);
[rel1,rel2,edge]=find(W_nnk_group);
new_Gr=graph(rel1,rel2,edge);
LWidths = 5*new_Gr.Edges.Weight/max(new_Gr.Edges.Weight);
figure(99),plot(new_Gr,'EdgeLabel',new_Gr.Edges.Weight,'LineWidth',LWidths)