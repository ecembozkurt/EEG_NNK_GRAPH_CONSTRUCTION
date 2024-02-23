function [alpha,beta,theta,gamma,delta]=subbands(eeg)
% delta :[0.5 4] hz
% theta : [4 8] hz
% alpha : [8 13] hz
% beta : [13 30] hz
% gamma : [30 45] hz
[h,w]=size(eeg);
alpha=zeros(h,w);
beta=zeros(h,w);
delta=zeros(h,w);
theta=zeros(h,w);
gamma=zeros(h,w);


for i=1:h
    S = eeg(i,:);
    waveletFunction = 'db8'; %OR 'sym8' OR 'coif5' OR 'db4';
    [C,L] = wavedec(S,8,waveletFunction);
    %%Calculation The Coificients Vectors
    cD1 = detcoef(C,L,1);                   %NOISY
    cD2 = detcoef(C,L,2);                   %NOISY
    cD3 = detcoef(C,L,3);                   %NOISY
    cD4 = detcoef(C,L,4);                   %NOISY
    cD5 = detcoef(C,L,5);                   %GAMA
    cD6 = detcoef(C,L,6);                   %BETA
    cD7 = detcoef(C,L,7);                   %ALPHA
    cD8 = detcoef(C,L,8);                   %THETA
    cA8 = appcoef(C,L,waveletFunction,8);   %DELTA
    
    %%%%Calculation the Details Vectors
    D1 = wrcoef('d',C,L,waveletFunction,1); %NOISY
    D2 = wrcoef('d',C,L,waveletFunction,2); %NOISY
    D3 = wrcoef('d',C,L,waveletFunction,3); %NOISY
    D4 = wrcoef('d',C,L,waveletFunction,4); %NOISY
    D5 = wrcoef('d',C,L,waveletFunction,5); %GAMMA
    D6 = wrcoef('d',C,L,waveletFunction,6); %BETA
    D7 = wrcoef('d',C,L,waveletFunction,7); %ALPHA
    D8 = wrcoef('d',C,L,waveletFunction,8); %THETA
    A8 = wrcoef('a',C,L,waveletFunction,8); %DELTA
    
    alpha(i,:)=D7;
    delta(i,:)=A8;
    theta(i,:)=D8;
    beta(i,:)=D6;
    gamma(i,:)=D5;
end
