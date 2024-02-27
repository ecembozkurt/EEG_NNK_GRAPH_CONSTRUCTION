function [W,error] = nnk_inverse_kernel_graph(G, mask, knn_param, reg)
% Constructs NNK graph
% G = similairty or kernel matrix
% mask - Neighbors to consider at each node (each row corresponds to neighbors)
% k_choice - Not used. For adaptive kernels
% reg - for removing weights smaller than reg value

if nargin < 4
    reg=1e-6;
end
n = sum(any(mask,2)); %min(size(G)); knn_mask, nonzero rows
neighbor_indices = zeros(min(size(G)), knn_param);
weight_values = zeros(min(size(G)),knn_param);
error_values = zeros(min(size(G)),knn_param);
nonzero_rows_only=[find(any(mask,2))]' ;
for i =nonzero_rows_only %nonzero rows only
    nodes_i = find(mask(i,:)); % only consider similar enough neighbors
    nodes_i(nodes_i==i) = [];
    G_i = full(G(nodes_i, nodes_i));% removing i-th row and column from G
    g_i = full(G(nodes_i,i)); % removing the i-th elem from
    C1=cond(G_i);
    if C1 >100
        C1
    end

    qp_output = nonnegative_qp_solver(G_i, g_i, reg, g_i);
    qpsol = qp_output.xopt;

    neighbor_indices(i,1:length(nodes_i)) = nodes_i; %changed
    weight_values(i,1:length(nodes_i)) = qpsol; %changed
    error_values(i,1:length(nodes_i)) = (G(i,i) - 2*qpsol'*g_i + qpsol'*G_i*qpsol); %changed

end

%%
%if C1<100
W=zeros(min(size(G)), min(size(G)));
error=zeros(min(size(G)), min(size(G)));
%     row_indices = repmat(nonzero_rows_only, n-1,1)';
%     nonzero_neighbor_indices= neighbor_indices(nonzero_rows_only,1:(n-1));
%     weight_values= weight_values(nonzero_rows_only,1:(n-1));
%     error_values=error_values(nonzero_rows_only,1:(n-1));

for iip=nonzero_rows_only
    nonzero_neighbor_num=sum(all(error_values(iip,:),1));
    W(iip,neighbor_indices(iip,1:nonzero_neighbor_num))=weight_values(iip,1:nonzero_neighbor_num);
    error(iip,neighbor_indices(iip,1:nonzero_neighbor_num))=error_values(iip,1:nonzero_neighbor_num);
end

%W = sparse(row_indices(:),  nonzero_neighbor_indices(:), weight_values(:), min(size(G)),min(size(G)));
%error = sparse(row_indices(:), nonzero_neighbor_indices(:), error_values(:), min(size(G)), min(size(G)));

W=sparse(W);
error=sparse(error);
%    W(error > error') = 0;
% W = max(W, W');
%else
%error=0;W=0;
% return
%end
%Ecem: track error


