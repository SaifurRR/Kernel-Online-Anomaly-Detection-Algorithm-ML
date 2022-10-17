%%% Experiment 14 %%%
%%% Investigate performance of KOAD as a function of thresholds nu1 and nu2%%%

clc; clear; close all;

% Load value of optimal sigma choice %
load C:\Users\tarem\LabComp\Workspace\IntruderDet\Indoor_IIUM\Exps\Exp2.mat X sigma_star T_train anomalies;

[T D] = size(X);

kernelChoice = 2;
sigma = sigma_star;
%%%%%%%%%%%%%%%%%%%%
sigma = 2*sigma_star;
%%%%%%%%%%%%%%%%%%%%

% nu1 = 0.10;
% nu2 = 0.60;

el = 1;  %Parameters for resolving orange alarm
epsilon = 0.2; 
L = 6;  %Parameters for dropping obsolete elements
d = 0.9;


index_actual=zeros(T,1); index_actual(anomalies) = 1;
index_actual([1:T_train]) = 0; %Ignore the training period
actual = sum(index_actual);

deltaStore_save={};

clear det_KOAD false_KOAD;
nu1_set = [0.01 0.05:0.05:0.95 0.99]
for n1=1:length(nu1_set)
    nu1 = nu1_set(n1)
    nu2_set = [nu1:0.05:0.99];
    for n2=1:length(nu2_set)
        nu2 = nu2_set(n2);

        [Red1 Red2 deltaStore Error] = KOAD(X, nu1, nu2, kernelChoice, sigma, d, L, epsilon, el);
        deltaStore_save{n1,n2} = deltaStore;
        index_det_KOAD = zeros(T,1); index_det_KOAD(Red1) = 1; index_det_KOAD(Red2) = 1; index_det_KOAD(1:T_train) = 0; %Count both Red1 and Red2 as detected
        det_KOAD(n1,n2) = length(find(index_actual==1 & index_det_KOAD==1));
        index_false_KOAD = zeros(T,1); index_false_KOAD(find(index_actual==0 & index_det_KOAD==1)) = 1;
        false_KOAD(n1,n2) = sum(index_false_KOAD);
    end %for n2=1:length(nu2_set)
    figure;
    plot(false_KOAD(n1,:)./(T-length(anomalies)), det_KOAD(n1,:)./length(anomalies));
    axis([0 1 0 1]);
    title(nu1);
end %for n1=1:length(nu1_set)

detRate_KOAD = det_KOAD./actual;
FDR_KOAD = false_KOAD./(det_KOAD+false_KOAD);
false_alarms_KOAD = false_KOAD./(T-actual);


% save C:\Users\tarem\LabComp\Workspace\IntruderDet\Indoor_IIUM\Exps\Exp14.mat;

n1=13; nu1=nu1_set(n1)
figure; clf;
plot(FDR_KOAD(n1,:), detRate_KOAD(n1,:), 'b');
xlabel('FDR'); ylabel('P_D');
title('KOAD');

figure; clf;
plot(false_alarms_KOAD(n1,:), detRate_KOAD(n1,:), 'k');
xlabel('P_{FA}'); ylabel('P_D');
title('KOAD');

% save C:\Users\tarem\LabComp\Workspace\IntruderDet\Indoor_IIUM\Exps\Exp14.mat;
