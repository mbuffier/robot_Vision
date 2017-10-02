%% **** Task 1.2: Characterising the Performance of the VO Algorithm ******* 

% for this task, I used NUM_INITIALISE = 10 and WINDOW_SIZE = 10
[rte_all1, ate_all1, rte_window1, ate_window1, stand_dev1, meanMagSD1, ~, covAfterBA] = visual_odometry(10, 10, 0, 1, 1, 0) ;

% ************ draw the results of evolutions ************
% RET and ATE over time
Y = [1 2:500]' ;

fig1 = figure ;
set(fig1, 'Position', [0 0 700 600]),
subplot(2,1,1), 
plot(Y, rte_all1, 'g');
hold on
plot( Y,rte_window1, 'b'); 
hold off
legend('RTE for all position','RTE for window'), title('RTE comparison');

subplot(2,1,2), 
plot(Y, ate_all1,'g') ;
hold on
plot(Y, ate_window1, 'b') ;
hold off ; 
legend('ATE for all position','ATE for window' ), title('ATE comparison') ;
hold off ;

saveas(fig1,'task12_ate_rte.jpg') ;

% evolution of the covariance along the axis after the loop closure
stand_dev1 = stand_dev1(1:379, :) ;
covAfterBA = covAfterBA(1:379,:) ;
Y3 = [1 2:379] ;

fig2 = figure ;
set(fig2, 'Position', [0 0 700 600]),
subplot(3,1,1), 
plot(Y3, stand_dev1(:,1), 'r');
hold on, 
plot(Y3, covAfterBA(:,1), 'b');
hold off;
legend('x axis SD when camera i is added', 'x axis SD after BA'), title('Evolution of the x axis standard deviation for the camera position after loop closure');

subplot(3,1,2), 
plot(Y3, stand_dev1(:,2), 'r');
hold on, 
plot(Y3, covAfterBA(:,2), 'b');
hold off ;
legend('y axis SD when camera i is added', 'y axis SD after BA'), title('Evolution of the y axis standard deviation for the camera position after loop closure');

subplot(3,1,3), 
plot(Y3, stand_dev1(:,3), 'r');
hold on, 
plot(Y3, covAfterBA(:,3), 'b');
hold off; 
legend('z axis SD when camera i is added', 'z axis SD after BA'), title('Evolution of the z axis standard deviation for the camera position after loop closure');

saveas(fig2,'task12_covAxis.jpg') ;

% evolution of the mean of the standard deviation magnitude 
fig3 = figure ;
set(fig3, 'Position', [0 0 700 600]),
plot(Y, meanMagSD1, 'g');
legend('mean of the SD magnitude from the start to time i'), title('Mean of the standard deviation magnitude for the camera position over time');

saveas(fig3,'task12_SD.jpg') ;

%% **** Task 1.3: Characterising the Performance of the VO Algorithm ******* 
% test the algorithm for a sliding window of 5, 10, 15 and 20 
[rte_all5, ate_all5, ~, ~, ~, ~, time5,~] = visual_odometry(5, 5, 0,0, 1, 1) ;
[rte_all10, ate_all10, ~, ~, ~, ~, time10,~] = visual_odometry(10, 10, 0,0, 1, 1) ;
[rte_all15, ate_all15, ~, ~, ~, ~, time15,~] = visual_odometry(15, 15, 0,0, 1, 1) ;
[rte_all20, ate_all20, ~, ~, ~, ~, time20,~] = visual_odometry(20, 20, 0,0, 1, 1) ;
[rte_all25, ate_all25, ~, ~, ~, ~, time25,~] = visual_odometry(25, 25, 0,0, 1, 1) ;

% plot the curve to compare visually
fig4 = figure; 
set(fig4, 'Position', [0,0, 1000, 1000]);

subplot(3,1,1), 
plot(Y, rte_all5, 'g');
hold on
plot(Y, rte_all10, 'k');
hold on
plot(Y, rte_all15, 'r');
hold on
plot(Y, rte_all20, 'b');
hold on
plot(Y, rte_all25, 'y');
hold off
legend('WS = 5','WS = 10', 'WS = 15', 'WS = 20', 'WS = 25'), title('RTE comparison over time for all window sizes');

subplot(3,1,2), 
plot(Y, ate_all5, 'g');
hold on
plot(Y, ate_all10, 'k');
hold on
plot(Y, ate_all15, 'r');
hold on
plot(Y, ate_all20, 'b');
hold on
plot(Y, ate_all25, 'y');
hold off
legend('WS = 5','WS = 10', 'WS = 15', 'WS = 20',  'WS = 25'), title('ATE comparison over time for all window sizes');

subplot(3,1,3), 
plot(Y, time5, 'g');
hold on
plot(Y, time10, 'k');
hold on
plot(Y, time15, 'r');
hold on
plot(Y, time20, 'b');
hold on
plot(Y, time25, 'y');
hold off
legend('WS = 5','WS = 10', 'WS = 15', 'WS = 20',  'WS = 25'), title('Computational time for all window sizes');

saveas(fig4,'task13.jpg') ;

% plot a comparison of all criteria in % for the different window sizes
ate_comparison = [ate_all5(end,1), ate_all10(end,1), ate_all15(end,1), ate_all20(end,1),ate_all25(end,1)] ;
ate_comparison = ate_comparison./max(ate_comparison(:)) ;

rte_comparison = [rte_all5(end,1), rte_all10(end,1), rte_all15(end,1), rte_all20(end,1),rte_all25(end,1)] ;
rte_comparison = rte_comparison./max(rte_comparison(:)) ;

times = [sum(time5,1), sum(time10,1), sum(time15,1), sum(time20,1),sum(time25,1)] ;
times = times./max(times(:)) ;

Y2 = [5:5:25] ;
fig42 = figure ;
subplot(3,1,1), 
plot(Y2,ate_comparison, '-*'), title('ATE result comparison (in %) in function of different window size'), 
xlabel('Window Size') ;

subplot(3,1,2), 
plot(Y2,rte_comparison, '-*'), title('RTE result comparison  (in %) in function of different window size'), 
xlabel('Window Size') , ylabel('% of each with respect to the maximum result');

subplot(3,1,3), plot(Y2,times, '-*'), title('Total time computation result (in %) in function of different window size'), 
xlabel('Window Size') ;

saveas(fig42,'task13_%.jpg') ;

%% ******************* Task 1.4: Effect of the GPS ****************** 

% comparaison by adding the GPS prior at frame 1, 250 and 500
[rte_all, ate_all, ~, ~, ~, ~, ~, ~] = visual_odometry(15, 15, 0, 0, 1, 0) ;
[rte_all_GPS, ate_all_GPS, ~, ~, ~, ~, ~, ~] = visual_odometry(15, 15, 1, 0, 1, 0) ;

fig5 = figure ; 
set(fig5, 'Position', [0,0, 900, 700]);

subplot(2,1,1), 
plot(Y, rte_all, 'g');
hold on
plot(Y, rte_all_GPS, 'b');
hold off
legend('without GPS','with GPS'), title('RTE comparaison over time with and without GPS');

subplot(2,1,2), 
plot(Y, ate_all, 'g');
hold on
plot(Y, ate_all_GPS, 'b');
hold off
legend('without GPS','with GPS'), title('ATE comparaison  over time with and without GPS');

saveas(fig5,'task14.jpg') ;

