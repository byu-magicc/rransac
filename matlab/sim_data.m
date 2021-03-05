se2_data = csvread("output2_se2.csv",0,1);
r2_data1 = csvread("output_r21.csv",0,1);
r2_data2 = csvread("output_r22.csv",0,1);
r2_data3 = csvread("output_r23.csv",0,1);
r2_data4 = csvread("output_r24.csv",0,1);
r2_data5 = csvread("output_r25.csv",0,1);
r2_data6 = csvread("output_r26.csv",0,1);
r2_data7 = csvread("output_r27.csv",0,1);
r2_data8 = csvread("output_r28.csv",0,1);
r2_data9 = csvread("output_r29.csv",0,1);
r2_data10 = csvread("output_r210.csv",0,1);
r2_data = zeros(14,101);

for i = 1:14
  r2_data(i,:) = (r2_data1(i,:) + r2_data2(i,:) + r2_data3(i,:) + r2_data4(i,:) + r2_data5(i,:) + r2_data6(i,:) + r2_data7(i,:) + r2_data8(i,:) + r2_data9(i,:) + r2_data10(i,:))/10;   
end
alpha = 0.5;
fig = figure(1),clf;
size = 30;
colormap(gray)
subplot(1,2,1);
plot(r2_data(5,:),r2_data(1,:),'*', 'color', [0,0,0] );
hold on
plot(se2_data(5,:),se2_data(1,:), 'color', [0,0,0]+alpha);
xlabel('time (s)');
ylabel('NVT');
legend('LTI-CV','SE2');
set(gca,'FontSize',size);

subplot(1,2,2);
plot(r2_data(5,:),r2_data(3,:),'*', 'color', [0,0,0]);
hold on
plot(se2_data(5,:),se2_data(3,:), 'color', [0,0,0]+alpha);
xlabel('time (s)');
ylabel('NMT');
legend('LTI-CV','SE2');
set(gca,'FontSize',size);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'results','-dpdf','-r0')



