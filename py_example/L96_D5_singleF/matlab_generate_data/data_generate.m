clear
T_total = 0:0.025:100;
data_initial = randn(5+1,1);%[0.80; 0.95; 0.71; 0.24; 0.63];
data_initial(end)=8.17;
[T,Y] = ode45(@lorenz_origin,T_total,data_initial);
save('data_true.mat','Y')
dlmwrite('clean_data.dat',Y(1000:2500,:),'delimiter',' ','precision','%e');
Y = Y + 0.5.*randn(size(Y));
save('data.mat','Y')
dlmwrite('twin_data.dat',Y(1000:2500,:),'delimiter',' ','precision','%e');