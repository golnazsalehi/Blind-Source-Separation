%% MOD
N0 = 3;
x1 = X;
[n,m] = size(D);
d = randn(n,m);
d = d./sqrt(sum(d.^2));
[n,m] = size(S);
MOD = zeros(1,100);
for ok =1:100    
    % source finder
    posOMP=zeros(1,N0);
    sOMP=zeros(n,m);
    for L = 1:m
        sOMP(:,L) = S_finder(X(:,L),x1(:,L),d,sOMP(:,L),posOMP,N0);
    end
	% dictionary finder
    d = x1*pinv(sOMP);
    d = d./sqrt(sum(d.^2));
    % Representation Error of MOD
    X_hat = d*sOMP;
    MOD(ok) = trace((X-X_hat)*(X-X_hat)')/trace(X*X');   
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
fprintf("\n Successful Recovery Rate of MOD: %.3f \n",SRR);
figure(1)
plot(MOD)
xlabel('Iterations')
xlabel('Representation Error of MOD')
%% K-SVD
N0 = 3;
x1 = X;
[n,m] = size(D);
d = randn(n,m); 
d = d./sqrt(sum(d.^2));
[N,M] = size(S);
K_SVD = zeros(1,100);
for ok =1:100    
    % source finder
    posOMP=zeros(1,N0);
    sOMP=zeros(N,M);
    for L = 1:M
        sOMP(:,L) = S_finder(X(:,L),x1(:,L),d,sOMP(:,L),posOMP,N0);
    end
	% dictionary finder
    for i = 1:50
        index = find(sOMP(i,:) ~= 0);
        
        num = 1:50;
        num(i)=[];
        X_r = X - d(:,num)*sOMP(num,:);
        X_m = X_r(:,index);
        
        [U,LANDA,V] = svd(X_m);
        
        [m,n] = size(LANDA);
        L = min(m,n);
        [landa,pos]=sort(diag(LANDA(1:L,1:L)),'descend');
        d(:,i) = U(:,pos(1));
        sOMP(i,index) = landa(1)* (V(:,pos(1)))';
    end
    % Representation Error of K_SVD
    X_hat = d*sOMP;
    K_SVD(ok) = trace((X-X_hat)*(X-X_hat)')/trace(X*X');   
end
figure(1)
hold on
plot(K_SVD)

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
fprintf("\n Successful Recovery Rate of K-SVD: %.3f \n",SRR);
%% functions
function sOMP = S_finder(x,x1,D,sOMP,posOMP,N0)
   for i=1:N0
        ro=x1'*D;
        [~,posOMP(i)]=max(abs(ro));
        if i>1
             Dsub=D(:,posOMP(1:i));
             sOMP(posOMP(1:i),:)=pinv(Dsub)*x;
             x1=x-D*sOMP;
        else
             sOMP(posOMP(1))=ro(posOMP(1));
            x1=x1-sOMP(posOMP(1))*D(:,posOMP(1)); 
        end
   end
end