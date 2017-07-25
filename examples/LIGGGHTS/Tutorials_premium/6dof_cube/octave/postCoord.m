#!/usr/bin/octave --silent
% clear workspace

clear
close all
more off

set(0,"defaultlinelinewidth",2);

% --- user input ---
simResultsPath    =  '../'; % '../archive/v1_withRestart/results'; %
simResultsPattern = 'coord.txt';
% ------------------

filename = fullfile(simResultsPath,simResultsPattern);

% load data
data = load(filename);
time = data(:,1);
pos = data(:,2:4);
quat = data(:,5:8);

% plots
figure(1)
subplot(2,1,1);
plot(time,pos);
ylabel('xcm of cube in m');
legend('x','y','z');

subplot(2,1,2);
plot(time,quat);
xlabel('Time in s');
ylabel('quat of cube');
legend('q0','q1','q2','q3');