% load('results_10_28(small hw).mat')
% figure(1)
% plot(Costz,"*")
% hold on 
% plot(Costskip,"*")
% plot(Costsplit,"*")
% plot(costprop, "*")
% plot(Costnorm,"*")
% plot(costnorm2,"*")
% legend('Zaid split', 'skip','skip+split', 'proposed', 'no control-4', 'no control-2')
% xlabel('Simulation time (hour)')
% ylabel('Average travel time')
% fprintf('split: mean - %f, std - %f \n', mean(Costsplit), std(Costsplit))
% fprintf('skip: mean - %f, std - %f \n', mean(Costskip), std(Costskip))
% fprintf('zaid: mean - %f, std - %f \n', mean(Costz), std(Costz))
% fprintf('norm-4: mean - %f, std - %f \n', mean(Costnorm), std(Costnorm))

%colA = ["#FF0000",'#FFFF00','#00EAFF','#AA00FF','#FF7F00','#BFFF00','#0095FF','#FF00AA','#FFD400','#6AFF00','#0040FF','#EDB9B9','#B9D7ED','#E7E9B9','#DCB9ED','#B9EDE0','#8F2323','#23628F','#8F6A23','#6B238F','#4F8F23','#000000','#737373','#CCCCCC'];
colA = ["#00FF00", "#FFA500", "#0000FF"];
close all
figure(1)
hold on
plot(t_n, P_veh_n, 'Color', colA(3), LineWidth=1.3)%, DisplayName='No-control', LineWidth=1.3)
plot(t_z, P_veh_z, 'Color', colA(2),LineWidth=1.3)%, DisplayName='Benchmark', LineWidth=1.3)
plot(t_p, P_veh_p, 'Color', colA(1), LineWidth=1.3)%, DisplayName='Proposed', LineWidth=1.3)


plot(t_z, P_wait_z, 'HandleVisibility', 'off', Color=colA(2), LineStyle='--', LineWidth=1.3)
plot(t_p, P_wait_p, 'HandleVisibility', 'off', Color=colA(1), LineStyle='--', LineWidth=1.3)
plot(t_n, P_wait_n, 'HandleVisibility', 'off', Color=colA(3), LineStyle='--', LineWidth=1.3)

plot(t_p, P_walk_p, 'HandleVisibility', 'off', Color=colA(1), LineStyle=':', LineWidth=1.3)
legend({'No-control', 'Split policy (Benchmark)', 'Proposed'}, FontSize=12)
legend show
xlabel('Time (s)', 'FontSize',12)
ylabel('Cumulative count', FontSize=12)
ax = gca; 
ax.XGrid = 'off';
ax.YGrid = 'on';
xlim([4500 7500])
set(gca,"FontSize",12)
yticks(0:50:500)

%legend()

