function [W,error,ERR_ecem] = omp_nnk_inverse_kernel_graph(G, mask, knn_param, reg, perform_OMP)
% MP and OMP version of NNK
% Computes a sparse non-negative inverse of Gram matrix G
if nargin < 4
    reg=1e-6;
end
if nargin < 5
    perform_OMP = true;
end
n = min(size(G));
neighbor_indices = zeros(n, knn_param);
weight_values = zeros(n, knn_param);
error_values = zeros(n,knn_param);
% W = sparse(n,n);
% error = sparse(n,n);
% rcond_array = zeros(n,1);
% coherence_array = zeros(n,1);
% no_of_selected_nodes = zeros(n,1);
ERROR_TOL =reg;
%%
% sigma = k_choice;
% G = exp(-D./(2*sigma));
%%
  ERR_ecem=[];
for i = 1:n
    nodes_i = find(mask(i,:));
    nodes_i(nodes_i==i) = [];
    %     if isempty(nodes_i)
    %         fprint('No neighbors found for node %d', i);
    %     end
    G_i = G(nodes_i, nodes_i);
    g_i = G(nodes_i,i);
    total_no_of_neighbors = length(nodes_i);
    to_be_processed = g_i;
    to_be_processed_nodes = 1:total_no_of_neighbors;

    %     options = optimoptions('quadprog','Display','off', 'Algorithm','trust-region-reflective');
    % First step of OMP
    [omp_selection, omp_selection_node] = max(to_be_processed);
    %     omp_selection_node = find(to_be_processed == omp_selection);

    selected_nodes = [to_be_processed_nodes(omp_selection_node)];
    solution = [omp_selection]; % weight for first step is one
    prev_solution = solution;

    % Remove selected nodes
    to_be_processed(omp_selection_node) = [];
    to_be_processed_nodes(omp_selection_node) = [];

    iter = 0;
    prev_error = 1;
  
    while ~isempty(to_be_processed)
        [omp_selection, omp_selection_node] = max((to_be_processed - G_i(to_be_processed_nodes,selected_nodes)*solution));
        if (omp_selection < reg)
            break
        end
        selected_nodes = [selected_nodes to_be_processed_nodes(omp_selection_node)];
        %% Weights for selected atoms
        if perform_OMP

            qp_output = nonnegative_qp_solver(G_i(selected_nodes, selected_nodes), g_i(selected_nodes), reg, g_i(selected_nodes));
            solution = qp_output.xopt;
            
        else
            solution = [solution; omp_selection]; %use projection as weight % g_i(omp_selection_node)
        end

        % Remove selected node from to be processed nodes
        to_be_processed(omp_selection_node) = [];
        to_be_processed_nodes(omp_selection_node) = [];
        iter = iter+1;

        % Break out of the loop if error below threshold
        %  G(i,i) = 1
        err = (G(i,i) - 2*solution'*g_i(selected_nodes) + solution'*G_i(selected_nodes, selected_nodes)*solution);
        ERR_ecem(i,iter)=err;
        if abs(prev_error-err)./(prev_error)<=1e-3
            solution = prev_solution;
            selected_nodes = selected_nodes(1:end-1);
            disp('adaptive k:'),length(solution)
            break;
        end
        prev_error = err;

        prev_solution = solution;
        if err<ERROR_TOL
            break;
        end
    end
    %     reg_modified = min(max(g_i), reg); % Done to avoid graph being disconnected
    %     solution(solution < reg_modified) = 0;
    solution(solution < reg) = 0;
    %     G_i = G_i(selected_nodes, selected_nodes);
    %     rcond_array(i,1) = rcond(full(G_i));
    %     UG_i = triu(G_i,1);
    %     coherence_array(i,1) = max(UG_i(:));
    %     no_of_selected_nodes(i,1) = length(find(solution>0));
    %     error(nodes_i,i) = prev_error;
    %     W(nodes_i(selected_nodes),i) = solution;
    temp = zeros(1, knn_param);
    neighbor_indices(i,:) = nodes_i;
    temp(selected_nodes) = solution;
    weight_values(i,:) = temp;
    error_values(i,:) = prev_error;
end
%figure,plot(ERR_ecem(1,:),'LineWidth',1),xlabel('k'),ylabel('Error')
%%
% Symmetrization: can be done better e.g. (see Ming Yuan 2010)
% % disp(full(L));
% fprintf("Average conditioning of NNLS: %0.4f \n",mean(rcond_array));
% fprintf("Average coherence of NNLS: %0.4f \n",mean(coherence_array));
% fprintf("Error Statistics of NNLS: mean=%0.4f,std=%0.4f \n",mean(error),std(error));
% figure();
% subplot(4,1,1); plot(rcond_array)
% subplot(4,1,2); plot(coherence_array)
% subplot(4,1,3); plot(error)
% subplot(4,1,4); plot(no_of_selected_nodes)
%%
row_indices = repmat(1:n, knn_param,1)';
W = sparse(row_indices(:), neighbor_indices(:), weight_values(:),n,n);
error = sparse(row_indices(:), neighbor_indices(:), error_values(:), n, n);
%W(error > error') = 0;
%W = max(W, W');
%W = 0.5.*(W+ W');

end
