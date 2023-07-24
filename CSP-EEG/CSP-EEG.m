%% Setting all means of all channels to Zero
load hw5.mat
for i = 1:60
    for j = 1:30
        TrainData_class1(j,:,i) = TrainData_class1(j,:,i) - mean(TrainData_class1(j,:,i));
    end
end

for i = 1:60
    for j = 1:30
        TrainData_class2(j,:,i) = TrainData_class2(j,:,i) - mean(TrainData_class2(j,:,i));
    end
end

for i = 1:40
    for j = 1:30
        TestData(j,:,i) = TestData(j,:,i) - mean(TestData(j,:,i));
    end
end
%% Question1 : 
R_moving = zeros(30,30);
for i =1:60
    R_moving = R_moving + TrainData_class1(:,:,i)*(TrainData_class1(:,:,i))';
end
R_moving  = R_moving/60;
%--------------------------------
R_thinking = zeros(30,30);
for i =1:60
    R_thinking = R_thinking + TrainData_class2(:,:,i)*(TrainData_class2(:,:,i))';
end
R_thinking = R_thinking/60;
%--------------------------------
[W_CSP,S] = eig(R_thinking,R_moving);
for i = 1:30
    W_CSP(:,i) = W_CSP(:,i)/norm(W_CSP(:,i));
end
[val,col] = max(max(S));
w_1 = W_CSP(:,col);
[val,col] = min(max(S));
w_30 = W_CSP(:,col);
new_W_CSP = [w_1,w_30];
%--------------------------------
trial_49_c1 = TrainData_class1(:,:,49);
trial_49_c2 = TrainData_class2(:,:,49);
%--------------------------------
filtered_trial_49_c1 = new_W_CSP'*trial_49_c1;
filtered_trial_49_c2 = new_W_CSP'*trial_49_c2;
%--------------------------------
variance_c1=var(filtered_trial_49_c1');
variance_c2=var(filtered_trial_49_c2');

figure(1)
subplot(2,1,1)
plot(filtered_trial_49_c1(1,:))
hold on;
plot(filtered_trial_49_c1(2,:),'r')
title("Filtered 49th trial for class 1")
txt1 =['Variances =' num2str(variance_c1)];
text(100,3,txt1)

subplot(2,1,2)
plot(filtered_trial_49_c2(1,:))
title("Filtered 49th trial for class 2")
hold on;
plot(filtered_trial_49_c2(2,:),'r')
txt1 =['Variances =' num2str(variance_c2)];
text(100,10,txt1)

%[variance_c1 ;variance_c2]
%% Question2 : 
abs_w1 = abs(w_1);
abs_w30 = abs(w_30);
figure;
plot(abs_w1)
hold on 
plot(abs_w30)
legend('W1 ','W30')
%% Question3 : 
wCSP = zeros(30,14);
eigVals = max(S);
vals = eigVals;
for i = 1:7
    [w,c] = max(vals);
    find_index = eigVals==w;
    index = find(find_index,1);
    wCSP(:,i) = W_CSP(:,index);
    vals(:,c) = [];
end

i=14;
while i>7
    [w,c] = min(vals);
    find_index = eigVals==w;
    index = find(find_index,1);
    wCSP(:,i) = W_CSP(:,index);
    vals(:,c) = [];
    i = i-1;
end
moving_filtered = zeros(14,256,60);
for i =1:60
    moving_filtered(:,:,i) = wCSP'*TrainData_class1(:,:,i);
    
end
moving_features = zeros(14,1,60);
for i =1:60
    moving_features(:,:,i) = var(moving_filtered(:,:,i)');
end



thinking_filtered = zeros(14,256,60);
for i =1:60
    thinking_filtered(:,:,i) = wCSP'*TrainData_class2(:,:,i);
    
end
thinking_features = zeros(14,1,60);
for i =1:60
    thinking_features(:,:,i) = var(thinking_filtered(:,:,i)');
end
%-------------------------
think_mean = zeros(14,1);
for i =1:14
    Sum = 0;
    for j=1:60
        Sum = Sum + thinking_features(i,1,j);
    end
    think_mean(i,1) = Sum/60;
end

moving_mean = zeros(14,1);
for i =1:14
    Sum = 0;
    for j=1:60
        Sum = Sum + moving_features(i,1,j);
    end
    moving_mean(i,1) = Sum/60;
end

A = (moving_mean-think_mean)*(moving_mean-think_mean)';
%-------------------------
T_covariance = zeros(14,14);
for i =1:60
    difference = thinking_features(:,1,i) - think_mean;
    T_covariance = T_covariance + difference*difference';
end
T_covariance = T_covariance/60;


M_covariance = zeros(14,14);
for i =1:60
    difference = moving_features(:,1,i) - moving_mean;
    M_covariance = M_covariance + difference*difference';
end
M_covariance = M_covariance/60;

B = T_covariance + M_covariance;
%-------------------------
[W,S] = eig(A,B);
[val,col] = max(max(S));
W_LDA = W(:,col)/norm(W(:,col));
%-------------------------
M_1 = W_LDA'*think_mean;

M_2 = W_LDA'*moving_mean;
 
threshold = (M_1 + M_2)/2;
%% Question4 : 

test_filtered = zeros(14,256,40);
for i =1:40
    test_filtered(:,:,i) = wCSP'*TestData(:,:,i);
    
end
test_features = zeros(14,1,40);
for i =1:40
    test_features(:,:,i) = var(test_filtered(:,:,i)');
end
projections = zeros(1,40);
for i=1:40
    projections(:,i) = W_LDA'*test_features(:,1,i);
end

labels_1 = projections >= threshold;
labels_2 = 2*(projections < threshold);
labels = labels_1 + labels_2;
%% Question5 : 
figure;
scatter(1:40,labels,'m','LineWidth',3)
hold on
scatter(1:40,TestLabel,'c','filled')
xlabel('The trials')
ylabel('Projections and True Labels')
legend('Projections','True Labels')