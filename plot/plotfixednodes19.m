fig=figure;
%W=neighbor_mask;
%W=squeeze(W_M(1,:,:));
for i=1:ch_number
    
    %nexttile
    %subplot(8,8,i),
    
    p4=plot3(eeg_xyz(:,1),eeg_xyz(:,2),eeg_xyz(:,3),'*') ;  
   % set(fig,'Visible','off')
    hold all
    for k=1:ch_number
        if W(i,k)==0 || W(i,k)==Inf
        else
            cc= abs(W(i,k));
            p4=plot3([eeg_xyz(i,1),eeg_xyz(k,1)],[eeg_xyz(i,2),eeg_xyz(k,2)],[eeg_xyz(i,3),eeg_xyz(k,3)],'Color',[0,0,0]);
            grid on
            
            
            %p3.LineColor='g';
            %p3.LineWidth=10*abs(W(i,k)-min(min(W)))/abs(max(max(W)));
            
            p4.LineWidth=abs(W(i,k))*5;
            
        end
        
    end
    % title(sprintf('Channel: %d ' ,i))
end
hold off;
%G.Edges.LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
%p3.LineWidth=G.Edges.LWidths;

set(fig,'Position',[1000 0 1200 1200])
% 
% figure, imagesc(W)
% title('neighbor mask')
% title('Neighbor Mask')
% xlabel('channels'), ylabel('channels')
