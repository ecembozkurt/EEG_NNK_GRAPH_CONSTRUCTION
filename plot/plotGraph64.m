function fig=plotGraph64(W,s,ch_number,varargin)
node_locations=load('eeg_xyz.mat');
eeg_xyz=node_locations.eeg_xyz;
if ch_number==14
    eeg_xyz=load('emotiv14_xyz.mat').emotiv14_xyz;
end
if ch_number==19
    eeg_xyz=load('./data/eeg_xyz_1020.mat').eeg_xyz1020;
end
fig=figure('visible', 'off');
if nargin<=2 
    ch_number=64;
end
tiledlayout(int8(sqrt(ch_number)),int8(sqrt(ch_number)), 'Padding', 'none', 'TileSpacing', 'compact');
for i=1:ch_number
    
    nexttile
    subplot(4,5,i),
    
    p4=plot3(eeg_xyz(:,1),eeg_xyz(:,2),eeg_xyz(:,3),'*') ;  
    hold all
    for k=1:ch_number
        if W(i,k)==0 || W(i,k)==Inf
        else
            cc= abs(W(i,k))/10;
            p4=plot3([eeg_xyz(i,1),eeg_xyz(k,1)],[eeg_xyz(i,2),eeg_xyz(k,2)],[eeg_xyz(i,3),eeg_xyz(k,3)],'Color',[0,0.2,cc]);
            grid on
            p4.LineWidth=abs(W(i,k))/10;
            
        end
        
    end
    % title(sprintf('Channel: %d ' ,i))
end

%G.Edges.LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
%p3.LineWidth=G.Edges.LWidths;
if nargin==1
    suptitle('EEG signals 64 channel');
else
    suptitle(s)
end
set(fig,'Position',[1000 0 1200 1200])
set(fig,'Visible','off')

