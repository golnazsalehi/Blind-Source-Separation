%% Question1 - defining signals and noise and observation martix
fs = 10^6; %Hz
duration = 0.001; %sec
time = linspace(0,duration,fs*duration);
fc = 150 * 10^6; %Hz
c = 3*10^8; %m/c -> speed of light
k = 2*pi*fc/c;

M=10; 
d=(0:M-1)';


f1 = 20000;
s1 = @(t) exp(1i*2*pi*f1*t);

f2 = 10000;
s2 = @(t) exp(1i*2*pi*f2*t);

a1=exp(-1j*k*d*sind(10));
a2=exp(-1j*k*d*sind(20));

noise=randn(M,1000);

observation = a1*s1(time) + a2*s2(time) + noise;
%% Question2 - applying Beamforming to find the angles of the sources
[U,S,V] = svd(observation); 
% first three columns of U are in the mixture function space and first three rows of V are in the source space

a_H = @(o1) exp(-1j*k*d*sind(o1));             

theta= 0.1:0.01:90;

L = length(theta);

new_U = [U(:,1), U(:,2)]; %we have  two sources!

Save_theta = zeros(1,L);

for i= 1: L
        O = a_H(theta(i))'*new_U; 
        Save_theta(i) = norm(O);
end
findpeaks(Save_theta)
xlabel("number of column in theta matrix")
ylabel(" corresponding norm")
%% Question3 - applying MUSIC to find the angles of the sources
new_U=U(:,3:end);

Save_theta = zeros(1,L);

for i= 1: L     
        O = a_H(theta(i))'*new_U; 
        Save_theta(i) =1/norm(O)^2;
end

plot(theta,Save_theta)
xlabel("number of column in theta matrix")
ylabel(" corresponding value due to MUSIC method")
%% Question4 - applying Beamforming to find the frequencies of the sources
S_T = @(f,t) exp(1i*2*pi*f*t);

f = 5000:1:25000;

L = length(f);

new_V=V(:,1:2); 

Save_f = zeros(1,L);

for i= 1: L
        O = S_T(f(i),time)*new_V;
        Save_f(i) = norm(O);
end

findpeaks(Save_f)
xlabel("number of column in f matrix")
ylabel(" corresponding norm")
%% Question5 - applying MUSIC to find the frequencies of the sources
new_V= V(:,3:end);

Save_f = zeros(1,L);
for i= 1: L
        O = S_T(f(i),time)*new_V;
        Save_f(i) =1/norm(O)^2;
end

plot(f,Save_f)
xlabel("number of column in f matrix")
ylabel(" corresponding value due to MUSIC method")