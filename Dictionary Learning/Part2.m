%% Dictionary generation
while(true)
    D = randn(3,6);
    D = D./repmat(sqrt(sum(D.^2)),3,1);
    ItsOk = Mutual_Coherence(D);
    if ItsOk==1
        break
    end    
end
%% Source generation
S = zeros(6,1000);
for i = 1:1000
    Source = randi([1 6]);
    x = -5+10*rand(1,1);
    S(Source,i) = x;
end
%% Noise generation
Noise = randn(3,1000)*0.1;
%% Observation
X = D*S + Noise;
%% Scatter Plot of the Observations
scatter3 (X(1,:),X(2,:),X(3,:))
%% MOD
d = randn(3,6);
d = d./repmat(sqrt(sum(d.^2)),3,1);
x1 = X;
for ok =1:600    
    % source finder
    sMP=zeros(6,1000);
    for L = 1:1000
        sMP(:,L) = S_finder(x1(:,L),d,sMP(:,L));
    end
	% dictionary finder
    d = x1*pinv(sMP);
    d = d./repmat(sqrt(sum(d.^2)),3,1); 
end  
Saving_D = D;
count =0;
for i = 1:size(d,2)
    Max = [];
    for  j = 1:size(Saving_D,2)
        Max = [Max,abs((d(:,i))' * Saving_D(:,j))];
    end
    if(max(Max) > 0.99)
        count = count + 1;
        [~,indx] = max(Max);
        Saving_D(:,indx) = [];
    end
end
SRR = count/size(D,2)*100;
fprintf("\n Successful Recovery Rate: %.3f \n",SRR);
%% 
d = [d(:,4),-d(:,6),-d(:,3),-d(:,1),-d(:,2),d(:,5)];
sMP=zeros(6,1000);
for L = 1:1000
        sMP(:,L) = S_finder(x1(:,L),d,sMP(:,L));
end
Error_MOD = trace((sMP-S)*(sMP-S)')/trace(S*S');
fprintf("\n Error of MOD: %.3f \n",Error_MOD);
%% K-SVD
d = randn(3,6);
d = d./repmat(sqrt(sum(d.^2)),3,1);
x1 = X;
for ok = 1:1500   % source finder
    sMP=zeros(6,1000);
    for L = 1:1000
        sMP(:,L) = S_finder(x1(:,L),d,sMP(:,L));
    end
	% dictionary finder
    for i = 1:6
        index = find(sMP(i,:) ~= 0);
        
        num = 1:6;
        num(i)=[];
        X_r = X - d(:,num)*sMP(num,:);
        X_m = X_r(:,index);
        
        [U,LANDA,V] = svd(X_m);
        if isempty(LANDA) == 0
            [m,n] = size(LANDA);
            L = min(m,n);
            [landa,pos]=sort(diag(LANDA(1:L,1:L)),'descend');
            d(:,i) = U(:,pos(1));
            sMP(i,index) = landa(1)* (V(:,pos(1)))';
        end
    end
end 
Saving_D = D;
count =0;
for i = 1:size(d,2)
    Max = [];
    for  j = 1:size(Saving_D,2)
        Max = [Max,abs((d(:,i))' * Saving_D(:,j))];
    end
    if(max(Max) > 0.99)
        count = count + 1;
        [~,indx] = max(Max);
        Saving_D(:,indx) = [];
    end
end
SRR = count/size(D,2)*100;
fprintf("\n Successful Recovery Rate: %.3f \n",SRR);
%% 
d = [-d(:,2) -d(:,1) d(:,6) -d(:,4) d(:,5) d(:,3)];
sMP=zeros(6,1000);
for L = 1:1000
        sMP(:,L) = S_finder(x1(:,L),d,sMP(:,L));
end
Error_K_SVD = trace((sMP-S)*(sMP-S)')/trace(S*S');
fprintf("\n Error of MOD: %.3f \n",Error_K_SVD);
%% Functions
function ItsOk = Mutual_Coherence(D)
    C = D'*D;
    S = size(C,1);
    Eye = eye(S);
    C = abs(C- C.*Eye);
    mutual_coherence = max(max(C));
    if mutual_coherence > 0.9
        ItsOk = 0;
    else
        ItsOk = 1;
    end
end

function [s] = S_finder(x1,D,sMP)
    ro=x1'*D;
    [~,posMP]=max(abs(ro));
    sMP(posMP)=ro(posMP);
    s = sMP;
end


