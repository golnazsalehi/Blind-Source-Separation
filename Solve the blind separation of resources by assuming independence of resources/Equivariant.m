clear
load('hw14.mat')
X = A*S + Noise;
without_noise = A*S;
%%
warning('off', 'all')
B =[0.321,0.532,0.533; 
    0.227,0.41,0.282;
    0.321,0.821,0.81];
B = B./sqrt(sum(B.^2,2));
Y = B*X;
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
E = zeros(1,1000);
for i =1:1000
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
    
 
    rond_B = [psi1*X' ;psi2*X'; psi3*X']/1001  - inv(B)';
    B = (eye(3) - 0.1*rond_B*B')*B;
    B = B./sqrt(sum(B.^2,2));
	Y = B*X;
    
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
    
    Shat(1,:)=Shat(1,:)*2;
    Shat(2,:)=Shat(2,:)*(-3);
    Shat(3,:)=Shat(3,:)*(-1.75);
    E(i) = (norm(S-Shat,'fro'))^2 / (norm(S,'fro'))^2;


end
%%
figure
plot(E);
title('Convergance Diagram');

figure
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


subplot(3,1,1)
plot(S(1,:),'LineWidth',1); hold on
plot(Shat(1,:)*2,'LineWidth',1);
legend('Real','Estimated')
title('First source')

subplot(3,1,2)
plot(S(2,:),'LineWidth',1); hold on
plot(Shat(2,:)*(-3),'LineWidth',1);
legend('Real','Estimated')
title('Second source')

subplot(3,1,3)
plot(S(3,:),'LineWidth',1); hold on
plot(Shat(3,:)*(-1.75),'LineWidth',1);
legend('Real','Estimated')
title('Third source')
%%
Shat(1,:)=Shat(1,:)*2;
Shat(2,:)=Shat(2,:)*(-3);
Shat(3,:)=Shat(3,:)*(-1.75);
E = (norm(S-Shat,'fro'))^2 / (norm(S,'fro'))^2
