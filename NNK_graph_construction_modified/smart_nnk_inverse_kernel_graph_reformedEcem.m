%general purpose incremental NNK 
function [W, per_pixel_timing] = smart_nnk_inverse_kernel_graph_reformedEcem(X, wsz, sigma_f, sigma_d, reg)
%% NNK graph for images
if nargin < 5
    reg=1e-12;
end
%%
shape = size(X);
ch=shape(1); T=shape(2);
N = ch*T;

%% Mask position values

D=DistEuclideanPiotrDollar(X,X);
%%
g_i = exp(-D./(sigma_d));
%% NNK specific initial values
to_remove_flag = ones(mask_size-1,1, 'logical');
mask_intensities_init = zeros(mask_size-1, 1);
zero_flag = zeros(mask_size-1, 1, 'logical');

sigma_ratio = -sigma_f/sigma_d ;
%%
neighbor_indices = ones(N, mask_size - 1);
weight_values = zeros(N, mask_size - 1);
%% NNK Graph construction 
% profile on;
for row=1:ch
    for col=1:T

        row_values = row+mask_positions(:,2);
        col_values = col+mask_positions(:,1);
        mask_intensities = mask_intensities_init;
        curr_remove_flag = to_remove_flag;
        g_curr = g_i;
        itr = 1; % itr iterates in the sorted index
        if (row<wsz+1)||(row>ch-wsz-1)||(col<wsz+1)||(col>T-wsz-1)  
            outside_image_index = (row_values <=0 | row_values>ch | col_values<=0 | col_values>T);
            row_values(outside_image_index) = []; col_values(outside_image_index) = [];
            g_curr(outside_image_index)=[];
            curr_remove_flag(outside_image_index) = 0;
        else
            outside_image_index = zero_flag;
        end
        
        p_values = X(sub2ind(size(X),row_values, col_values));
        p_values = p_values(:);
        p_i = X(row,col); 
        
        mask_intensities(~outside_image_index) = p_values - p_i;
        %% Filter 
        while (itr < mask_size)
            curr_index = index(itr);
            curr_nodes_on_line = nodes_on_line_indexes{itr}; thresh = sigma_ratio.*thresh_multiplier{itr};
            % If a node had been removed by earlier node or doesn't have any 
            % node that it can prune, skip processing.
            if ~curr_remove_flag(curr_index) || isempty(curr_nodes_on_line)  
                itr = itr + 1;
                continue;
            end

            curr_intensity_diff = mask_intensities(curr_index);
            if curr_intensity_diff == 0
                to_remove = curr_nodes_on_line;
            else
                to_remove = curr_nodes_on_line((curr_intensity_diff - ...
                        mask_intensities(index(curr_nodes_on_line)))*curr_intensity_diff <= thresh);
            end
            curr_remove_flag(index(to_remove)) = 0;
            itr = itr + 1;
        end
        %%
        curr_remove_flag(outside_image_index) = [];
        mask_intensities(outside_image_index) = [];
        nbr_index = (row_values -1)*T + col_values;
        %%
        weights = (g_curr.*curr_remove_flag).*exp(-(mask_intensities.*curr_remove_flag).^2/(sigma_f));
        %%
        neighbor_length = length(nbr_index);
        node_index = (row-1)*T + col;
        neighbor_indices(node_index, 1:neighbor_length) = nbr_index;
        weight_values(node_index, 1:neighbor_length) = weights;
        %%
    end
end
% profile off; profile viewer;
row_indices = repmat(1:N, mask_size - 1,1)';
W = sparse(row_indices(:), neighbor_indices(:), weight_values(:), N, N);
W = max(W, W');
end


