clear
load('hw14.mat')
X = A*S + Noise;
[U,D] = eig(X*X');
W = ((D^(-1/2))*U');
Z = W*X;
%%
B = [1,0.9063,0.2607
    0.5841,0.8797,0.5944
    0.1078,0.8178,1];


kurty1 = zeros(1,100);
kurty2 = zeros(1,100);
kurty3 = zeros(1,100);


for j =1:100
    for i=1:3
        b = (B(i,:))';    
        save_b = b;
        y=b'*Z;
        b = mean(([y.^3;y.^3;y.^3].*Z),2)-3*b;
        kurty = mean(y.^4)-3*(mean(y.^2))^2;
        if i==2
            b1 = B(1,:)';
            b = (eye(3) - b1*b1')*b;
        elseif i==3
            b1 = B(1,:)';
            b2 = B(2,:)';
            b = (eye (3) - [b1 b2]*[b1 b2]')*b; 
        end
        b = b./sqrt(sum(b.^2,1));

        if i==1
            B = [b';B(2,:);B(3,:)];
            kurty1(j) = kurty;
            
        elseif i==2
            B = [B(1,:);b';B(3,:)];
            kurty2(j) = kurty;
        elseif i==3
            B = [B(1,:);B(2,:);b'];
            kurty3(j) = kurty;
        end     
    end
end
%%
BW = B*W;
BWA=BW*A
Shat =BW*X;

Shatd=Shat; Sd=S;
[~,r1]=max(abs(Shatd(1,:)*Sd'));
Sd(r1,:) =0;
Shat(r1,:)=Shatd(1,:);

[~,r2]=max(abs(Shatd(2,:)*Sd'));
Sd(r2,:) =0;
Shat(r2,:)=Shatd(2,:);

[~,r3]=max(abs(Shatd(3,:)*Sd'));
Sd(r3,:) =0;
Shat(r3,:)=Shatd(3,:);
%%
figure;
plot([abs(kurty1);abs(kurty2);abs(kurty3)]')
legend('kurt y1','kurt y2','kurt y3');

figure
subplot(3,1,1)
plot(S(1,:),'LineWidth',1); hold on
plot(Shat(1,:),'LineWidth',1);
legend('Real','Estimated')
title('First source')

subplot(3,1,2)
plot(S(2,:),'LineWidth',1); hold on
plot(-Shat(2,:),'LineWidth',1);
legend('Real','Estimated')
title('Second source')

subplot(3,1,3)
plot(S(3,:),'LineWidth',1); hold on
plot(-Shat(3,:),'LineWidth',1);
legend('Real','Estimated')
title('Third source')
%%
Shat(2,:)=-Shat(2,:);
Shat(3,:)=-Shat(3,:);
E = (norm(S-Shat,'fro'))^2 / (norm(S,'fro'))^2