clear
T_total = 0:0.05:1000;
data_initial = randn(7,1);%[0.80; 0.95; 0.71; 0.24; 0.63];
data_initial(4:7)=[5.0, 0.0797, 0.6898, 6.2723];
[T,Y] = ode45(@func_origin,T_total,data_initial);
save('data_true.mat','Y')
dlmwrite('clean_data.dat',Y(1000:2500,:),'delimiter',' ','precision','%e');
Y = Y + 0.5.*randn(size(Y));
save('data.mat','Y')
dlmwrite('twin_data.dat',Y(1000:2500,:),'delimiter',' ','precision','%e');
