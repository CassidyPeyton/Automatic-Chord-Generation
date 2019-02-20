A =   [0.5000,0.5000; 0.5000,0.5000];
B =   [0.4000,0.1000,0.5000; 0.1000,0.5000,0.4000];
Pi =  [0.5000,0.5000]';
O =   [3,1,1,3,2,3,2,2,2,3,2,2,2,2,2,3,3,1,1,2]';
[a,b,c] = BaumWelch_n(A,B,Pi,O,20);
th = 0.0000005;
[Anew,Bnew,Pinew,errors,count] = BaumWelch_th(A,B,Pi,O,th);
figure
plot(1:count,errors)
xlabel('iterations')
ylabel('Eu-divergence')
grid on

