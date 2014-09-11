
%%
clear;
IC = [2.8694    1.6217   35.4603];
T = [0:0.04:1000]';
[t, y ]= ode45(@Lorenz63, T, IC);
D = 0.25.*y(:,1);
save stimulus.dat D -ascii
xdata = [T D];
save current.dat xdata -ascii
T = [0:0.04:1000]';
current=load('current.dat');
P = [1.0 120 115 20 -12 0.3 10.6 10 30 1 5 25 15 0.1 0.4 5 -15 1 7];
IC= [4.2176    0.2393    0.0000    0.1562];
[t2,y2] = ode45(@HH3, T, IC, [], current, P);
Y = zeros(length(T),23);
Y(:,1:4) = y2;
Y(:,5:end) = repmat(P,length(T),1);
save true_data.dat Y -ascii
Y(:,1) = Y(:,1) + 0.5.*randn(size(Y(:,1)));
save twin_data.dat Y -ascii
%%
% figure(1); clf;
% subplot(2,1,1);
% plot(T,D);
% subplot(2,1,2);
% plot(T, y2(:,1));
% figure(2);clf;
% for i = 2:4
%     subplot(3,1,i-1);
%     plot(T, y2(:,i));
% end

% xdata = [T y2 D];
% zdata = [T y2(:,1) D];
% save HHtanhLong.txt zdata -ascii
% save xLong.dat xdata -ascii
