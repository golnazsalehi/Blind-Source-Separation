clear
load('hw14.mat')
X = A*S + Noise;
without_noise = A*S;
%%
[U,D] = eig(X*X');
W = ((D^(-1/2))*U');
Z = W*X;
B =[0.321,0.532,0.533; 
    0.227,0.41,0.282;
    0.321,0.821,0.81];

B = B./sqrt(sum(B.^2,2));
Y = B*Z;

mu = 0.01;
k = @(Y) [ones(1,1001)
     Y
     Y.^2
     Y.^3
     Y.^4
     Y.^5];
 k_prime = @(Y) [zeros(1,1001)
                 ones(1,1001)
                 2*Y
                 3*Y.^2
                 4*Y.^3
                 5*Y.^4];
Err = zeros(1,1000);             
for j=1:1000  
    k1 = k(Y(1,:));
    k2 = k(Y(2,:));
    k3 = k(Y(3,:));
    
    KK1 = (k1*k1')/1001; 
    KK2 = (k2*k2')/1001; 
    KK3 = (k3*k3')/1001; 
    
    k1_prime = mean(k_prime(Y(1,:)),2);
    k2_prime = mean(k_prime(Y(2,:)),2);
    k3_prime = mean(k_prime(Y(3,:)),2);
    
    theta1 = KK1\k1_prime;
    theta2 = KK2\k2_prime;
    theta3 = KK3\k3_prime;

    psi1 = theta1'*k1;
    psi2 = theta2'*k2;
    psi3 = theta3'*k3;
    
    df1 = psi1*Z';
    df2 = psi2*Z';
    df3 = psi3*Z';
 
    b1 = B(1,:)';
    b2 = B(2,:)';
    b3 = B(3,:)';
            
    b1 = b1 - mu* df1';
    b1 = b1./sqrt(sum(b1.^2,1));

    b2 = b2 - mu* df2';
	b2 = (eye(3) - b1*b1')*b2;
    b2 = b2./sqrt(sum(b2.^2,1));

    
    b3 = b3 - mu* df3';
    b3 = (eye (3) - [b1 b2]*[b1 b2]')*b3; 
    b3 = b3./sqrt(sum(b3.^2,1));
   
    B = [b1';b2';b3'];
    %----------------------------
    Err(j) = norm([df1;df2;df3])/1001; 
    Y = B*Z;
    
end   
%%
figure
plot(Err);
title('Convergance Diagram');


Shat =B*X;

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

figure
subplot(3,1,1)
plot(S(1,:),'LineWidth',1); hold on
plot(-(2/3)*Shat(1,:),'LineWidth',1);
legend('Real','Estimated')
title('First source')

subplot(3,1,2)
plot(S(2,:),'LineWidth',1); hold on
plot(Shat(2,:),'LineWidth',1);
legend('Real','Estimated')
title('Second source')

subplot(3,1,3)
plot(S(3,:),'LineWidth',1); hold on
plot(-Shat(3,:),'LineWidth',1);
legend('Real','Estimated')
title('Third source')
%%

Shat(3,:)=-Shat(3,:);
Shat(1,:)=-(2/3)*Shat(1,:);
E = (norm(S-Shat,'fro'))^2 / (norm(S,'fro'))^2