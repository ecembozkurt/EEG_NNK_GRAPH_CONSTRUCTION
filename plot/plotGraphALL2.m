function fig=plotGraphALL2(W,s,ch_number,varargin)
if nargin>3
    mDelay=varargin{1};
else 
    mDelay=0;
end

if ch_number==19
    eeg_xyz=load('./data/eeg_xyz_1020.mat').eeg_xyz1020;
end
fig=figure;%('visible', 'off');
%set(fig,'Visible','off')
set(0, 'currentfigure', fig);
if nargin<=2
    ch_number=64;
end

p4=plot3(eeg_xyz(:,1),eeg_xyz(:,2),eeg_xyz(:,3),'*') ;

eeg_xyz=repmat(eeg_xyz,(mDelay+1),1);
for i=1:ch_number
    if i<=ch_number
    bbb = num2str(i); ccc = cellstr(bbb);
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
    text(eeg_xyz(i,1)+dx, eeg_xyz(i,2)+dy, eeg_xyz(i,3),ccc);
    hold all
    end
    hold all
    for k=1:(mDelay+1)*ch_number
        if W(i,k)<10^-3 || W(i,k)==Inf
        else
            cc= abs(W(i,k));
            p4=plot3([eeg_xyz(i,1),eeg_xyz(k,1)],[eeg_xyz(i,2),eeg_xyz(k,2)],[eeg_xyz(i,3),eeg_xyz(k,3)],'Color',[0.3*fix((k-1)/ch_number +1),1-0.3*fix((k-1)/ch_number +1),1]);
            
            grid on
            p4.LineWidth=abs(W(i,k))*4;
        end
    end
end
hold off;
set(fig,'Position',[1000 0 1200 1200])
title(s)


