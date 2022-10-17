clc
clear all
% tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p=xlsread('ICU_Patient_v4.xls');
% X=double(xlsread('ICU_Patient_v4.xls'));
% % Expt was performed with v4 excel file deltastore_medi

% p=xlsread('ICU_Patient_v5.xls');
% X=double(xlsread('ICU_Patient_v5.xls'));
% deltastore_medi_ions

p=xlsread('ICU_Patient_realdata_v1.xls');
X=double(xlsread('ICU_Patient_realdata_v1.xls'));
% deltastore_medi_pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% p=xlsread('ICU_Patient1.xls');
% X=double(xlsread('ICU_Patient1.xls'));
% p=xlsread('PATIENT TABLE_PR.xls');
% X=double(xlsread('PATIENT TABLE_PR.xls'));
%p=xlsread('PATIENT_TABLE_SBP.xlsx');
%X=double(xlsread('PATIENT_TABLE_SBP.xlsx'));
%p=xlsread('PATIENT_TABLE_DBP.xlsx');
%X=double(xlsread('PATIENT_TABLE_DBP.xlsx'));
%p=xlsread('PATIENT_TABLE_MBP.xlsx');
%X=double(xlsread('PATIENT_TABLE_CVP.xlsx'));
%p=xlsread('PATIENT_TABLE_SBP.xlsx');
%X=double(xlsread('PATIENT_TABLE_CVP.xlsx'));
%p=xlsread('PATIENT_TABLE_CO.xlsx');
%X=double(xlsread('PATIENT_TABLE_CO.xlsx'));
%p=xlsread('PATIENT_TABLE_CI.xlsx');
%X=double(xlsread('PATIENT_TABLE_CI.xlsx'));
%p=xlsread('PATIENT_TABLE_PASP.xlsx');
%X=double(xlsread('PATIENT_TABLE_PASP.xlsx'));
%p=xlsread('PATIENT_TABLE_PADP.xlsx');
%X=double(xlsread('PATIENT_TABLE_PADP.xlsx'));
%p=xlsread('PATIENT_TABLE_MPAP.xlsx');
%X=double(xlsread('PATIENT_TABLE_MPAP.xlsx'));
%p=xlsread('PATIENT_TABLE_PCWP.xlsx');
%X=double(xlsread('PATIENT_TABLE_PCWP.xlsx'));
%p=xlsread('PATIENT_TABLE_LAP.xlsx');
%X=double(xlsread('PATIENT_TABLE_LAP.xlsx'));
%p=xlsread('PATIENT_TABLE_LVEDP.xlsx');
%X=double(xlsread('PATIENT_TABLE_LVEDP.xlsx'));
%p=xlsread('PATIENT_TABLE_DELTAT.xlsx');
%X=double(xlsread('PATIENT_TABLE_DELTAT.xlsx'));
%p=xlsread('PATIENT_TABLE_MEANST.xlsx');
%X=double(xlsread('PATIENT_TABLE_MEANST.xlsx'));
%p=xlsread('PATIENT_TABLE_MPAP.xlsx');
%X=double(xlsread('PATIENT_TABLE_MPAP.xlsx'));
%p=xlsread('PATIENT_TABLE_PADP.xlsx');
%X=double(xlsread('PATIENT_TABLE_PADP.xlsx'));
%p=xlsread('PATIENT_TABLE_PASP.xlsx');
%X=double(xlsread('PATIENT_TABLE_PASP.xlsx'));
%p=xlsread('PATIENT_TABLE_PCWP.xlsx');
%X=double(xlsread('PATIENT_TABLE_PCWP.xlsx'));
%p=xlsread('PATIENT_TABLE_TEMPC.xlsx');
%X=double(xlsread('PATIENT_TABLE_TEMPC.xlsx'));
%p=xlsread('PATIENT_TABLE_TEMPP.xlsx');
%X=double(xlsread('PATIENT_TABLE_TEMPP.xlsx'));
%p=xlsread('Normal Patient _ ETCO2.xls');
%X=double(xlsread('Normal Patient _ ETCO2.xls'));
%p=xlsread('Normal Patient _ SPO2.xls');
%X=double(xlsread('Normal Patient _ SPO2.xls'));
%p=xlsread('Normal Patient_PaCO2.xls');
%X=double(xlsread('Normal Patient_PaCO2.xls'));
%p=xlsread('Normal Patient_SCV02 SVO2.xls');
%X=double(xlsread('Normal Patient_SCV02 SVO2.xls'));
%p=xlsread('Normal Patient_SVRI.xls');
%X=double(xlsread('Normal Patient_SVRI.xls'));

[T P]=size(X);


trueint=zeros(1,T);
act=[17 21 29 37 64 72 77 84 85 88 89];
trueint(act)=1;

%if nargin < 11 = 1; end %Needed only if using Gaussian kernel
gamma = 1; %Forgetting factor
r = 1;  %Parameters for resetting P
R = 10000; 
el = 10;  %Parameters for resolving orange alarm
epsilon = 0.2; 
L = 100;  %Parameters for dropping obsolete elements
d = 0.9; 
kernelChoice = 2;
sigma = 0.03;

X = X./repmat(sqrt(sum(X.*X,2)+eps),1,size(X,2)); %normalize to unit circle (i.e. divide by norm)
Y = sum(X,2); %Add after normalizing
[T P] = size(X);
sumdec=0;

nu1 = 0.05
     for nu2=0.01:0.01:1
%          nu2
        sumdec=sumdec+1;
         flagint=zeros(1,T);
%         reply = input('Do you want more? Y/N [Y]: ','s');
% if isempty(reply)
%     reply = 'Y';
%     close all
%  nu1 = 0.15 ; nu2 = 0.60; %Threshholds
%d = 0.9; L = 100; %Parameters for dropping obsolete elements
%epsilon = 0.20; el = 20; %Parameters for resolving orange alarm
%R = 10000; r = 1; %Parameters for resetting P
%gamma = 1; %Forgetting factor
%sigma = 1; %Needed only if using Gaussian kernel

Red1 = []; Red2 = [];%Clear alarms
Orange = []; x_Orange = []; %Store x in timesteps when Orange alarm is raised

% Initialize %
t = 1;
x = X(t,:)';
y = Y(t);
k11 = kernel(x, x,kernelChoice, sigma);
K_tilde = [k11];
K_tilde_inv = [1/k11];
Dictionary = [x];
index_m = [t]; %Keeps track of timesteps when elements are added (+2), deleted (-1) or no change to D (0); for debugging only
Orange = [Orange t];
x_Orange = [x_Orange x];
drop_index = [0];
P=[1]; %P=inv(A'A)
m=1;
m_t(t) = m; %Keep track of m, for debugging only
index_m(t) = 2; %index_m(t)=2 implies x(t) is being added to Dictionary
alpha = y(t)/k11;
deltaStore(1) = nu1+eps;  %For debugging

% % Evaluate y_hat %
y_hat = zeros(T,1);
Error = zeros(1,T);
for j=1:m
    y_hat(t) = y_hat(t) + alpha(j)*kernel(Dictionary(:,j),x,kernelChoice, sigma);
end %for j=1:m
Error(t) = (Y(t)-y_hat(t))/Y(t)*100;

clear Lambda lambda dotProd;
Lambda = kernel(Dictionary(:,j),x,kernelChoice, sigma);

%Keep track of all dot product (kernel) values; for debugging only
for j=1:m
    dotProd(t,j) = kernel(Dictionary(:,j),x,kernelChoice, sigma);
end %for j=1:m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:T
    %t=t+1;
    x = X(t,:)';
    y = Y(t);

    % Evaluate current measurement %
    k_tilde = zeros(m,1);
    for j=1:m
        k_tilde(j) = kernel(Dictionary(:,j),x,kernelChoice, sigma); %Computing k_tilde{t-1}
    end %for j=1:m
    a = K_tilde_inv*k_tilde;
    delta = kernel(x,x,kernelChoice, sigma) - k_tilde'*a;
%    deltaCheck = a'*K_tilde*a - 2*a'*k_tilde + kernel(x,x); %Verify delta; not part of algorithm
    deltaStore(t) = delta;  %Keep track of delta, for debugging only
    if t>L
        Lambda = [Lambda(2:end,1:end) ; ceil(k_tilde'-repmat(d,1,m))]; %Append with 1 or 0
    else
        Lambda = [Lambda ; ceil(k_tilde'-repmat(d,1,m))]; %Append with 1 or 0
    end %if t>L

    if (delta>=nu1 & delta<nu2) %Orange alarm, add to Dictionary
        x_Orange = [x_Orange x];
        Orange = [Orange t];
        Dictionary = [Dictionary x];
        drop_index = [drop_index 0];
        a_tilde = a;
        K_tilde_inv = [ (delta*K_tilde_inv+a_tilde*a_tilde') (-1*a_tilde) ; (-1*a_tilde') (1) ] / delta;
        K_tilde = [ K_tilde k_tilde ; k_tilde' kernel(x,x,kernelChoice, sigma) ];
        if t>L
            lambda = [zeros(L-1,1) ; 1];
        else
            lambda = [zeros(t-1,1) ; 1];
        end %if t>L
        Lambda = [Lambda lambda];

        a = [zeros(m-1,1) ; 1];
        P = [ P zeros(m,1) ; zeros(m,1)' gamma ]/gamma;
        alpha = [ (gamma^(-0.5)*alpha - a_tilde*(y-gamma^(-0.5)*k_tilde'*alpha)/delta)  ;  ((y-gamma^(-0.5)*k_tilde'*alpha)/delta) ];
        m=m+1;
        m_t(t) = m;
        index_m(t) = 2; %Element added to D in this timestep
    else %delta<nu1 or delta>=nu2, Dictionary unchanged
        if delta>nu2 %Red1 alarm
            Red1 = [Red1 t];
        end %if delta>nu2
        %K_tilde = K_tilde;
        %K_tilde_inv = K_tilde_inv;
        q = (P*a) / (gamma+a'*P*a);
        P = (1/gamma)*[P - q*a'*P];
        alpha = alpha + K_tilde_inv*q*(y-k_tilde'*alpha);
        %m = m;
        m_t(t) = m;
        index_m(t) = 0; %No change to D in this timestep
    end %delta > nu1 & delta < nu2 


    %Keep track of all dot product (kernel) values; for debugging only
    for j=1:m
        dotProd(t,j) = kernel(Dictionary(:,j),x,kernelChoice, sigma);
    end %for j=1:m
    
    % Process previous orange alarm %
    if t>el & sum(Orange==t-el)==1 %means orange alarm at timestep t-el
        %Identify Dictionary element j corr. to the orange alarm at timestep t-el
        for j=1:m
            if x_Orange(:,Orange==t-el)==Dictionary(:,j)
                break;
            end %if x_Orange(:,Orange==t-el)==Dictionary(:,j)
        end %for j=1:m

        if sum(Lambda(end-el+1:end,j)) <= epsilon*el
            %Orange turns Red
            Red2 = [Red2  Orange(Orange==t-el)]; %Red2 alarm
            x_Orange(:,Orange==t-el) = [];
            Orange(Orange==t-el) = [];
            drop_index = [zeros(1,j-1) 1 zeros(1,m-j) ];
        else
            %Orange turns green
            x_Orange(:,Orange==t-el) = [];
            Orange(Orange==t-el) = [];
        end %if size(find(Lambda(end-el+1:end,j)<d),1) >= 0.80*25
    end %if t>el & sum(Orange==t-el)==1


    % Remove obsolete elements %
    for j=1:m
        %Dropping condition: kernel exists for past L timesteps, and is always < d
        if  ( t>L & sum(Lambda(1:end,j))==0 )
            drop_index(j) = 1;
        end %if  ( t>L & gt(Lambda(:,j),0) & lt(Lambda(:,j),d) )
    end %for j=1:m

    % DropElement(p) %
    if ( find(drop_index==1) & m>1 & t>r )
        t;
        p = min(find(drop_index==1)); %Drop Dictionary element # p
        %Reorganize K_tilde_p and K_tilde_inv_p, with p'th row/col moved to the end
        K_tilde = [ K_tilde(1:p-1,1:p-1) K_tilde(1:p-1,p+1:m) K_tilde(1:p-1,p)  ;  K_tilde(p+1:m,1:p-1) K_tilde(p+1:m,p+1:m) K_tilde(p+1:m,p)  ;  K_tilde(p,1:p-1) K_tilde(p,p+1:m) K_tilde(p,p) ];
        K_tilde_inv = [ K_tilde_inv(1:p-1,1:p-1) K_tilde_inv(1:p-1,p+1:m) K_tilde_inv(1:p-1,p)  ;  K_tilde_inv(p+1:m,1:p-1) K_tilde_inv(p+1:m,p+1:m) K_tilde_inv(p+1:m,p)  ;  K_tilde_inv(p,1:p-1) K_tilde_inv(p,p+1:m) K_tilde_inv(p,p) ];
        delta_p = 1/(K_tilde_inv(m,m));
        a_tilde_p = -delta_p*[K_tilde_inv(1:m-1,m)];
        K_tilde_inv = K_tilde_inv(1:m-1,1:m-1)-a_tilde_p*a_tilde_p'/delta_p;
        alpha = alpha - (1/delta_p)*[a_tilde_p*a_tilde_p'  -a_tilde_p  ;  -a_tilde_p'  1] *K_tilde*alpha;
        alpha = alpha(1:m-1);
        K_tilde = K_tilde(1:m-1,1:m-1);
        Dictionary(:,p) = [];
        drop_index(p) = [];
        Lambda(:,p) = [];
        dotProd(:,p) = []; %%%
        m=m-1;
        m_t(t) = m;
        index_m(t) = -1; %Element deleted from D in this timestep

        % Reset P %
        P = R*eye(m);
        for i_r=1:r
            k_tilde = zeros(m,1);
            for j=1:m
                k_tilde(j) = kernel(Dictionary(:,j),X(t-i_r,:)',kernelChoice, sigma); %Computing k_tilde{t-1}
            end %for j=1:m
            a = K_tilde_inv*k_tilde;
            q = (P*a) / (gamma+a'*P*a);
            P = (1/gamma)*[P - q*a'*P];
            alpha = alpha + K_tilde_inv*q*(Y(t-i_r)-k_tilde'*alpha);
        end %for i_r=1:r-1
     end %if ( find(drop_index==1) & m>1 & t>r )

     % Evaluate y_hat %
     for j=1:m
         y_hat(t) = y_hat(t) + alpha(j)*kernel(Dictionary(:,j),x);
     end %for i=1:m
     Error(t) = (Y(t)-y_hat(t))/Y(t)*100;

end %for t=2:T


    Red1_out = Red1; Red2_out = Red2; 

    Red1_out = Red1; Red2_out = Red2; deltaStore_out = deltaStore;

    Red1_out = Red1; Red2_out = Red2; deltaStore_out = deltaStore; Error_out = Error;


 
       

Red1
Red2
flagint(Red1)=1;
flagint(Red2)=1;
flagint(1)=0;
flagint;
detected=bitand(flagint,trueint);
false=bitxor(flagint,trueint);
falsem=false;
false(act)=0;
missed=bitxor(falsem,false);
mis(sumdec)=sum(missed);
dec(sumdec)=(sum(detected)/11)*100;
fal(sumdec)=(sum(false)/(100-11))*100;
figure(1)
 scatter(sort(fal),sort(dec));
% scatter(fal,dec);
% plot(sort(fal),sort(dec));

    end
%      end
figure(2)
% stem(deltaStore_out)
% xlabel('timesteps')
% ylabel('deltastore_out')


g = stem(deltaStore_out, 'k');
set(g, 'LineWidth', 1);
hold on;
g = stem(act,deltaStore_out(act), 'r', 'filled');
xlabel('timesteps')
ylabel('deltastore_out')
% th=[nu1;nu2];
% legend(th)
% legend('nu1= ', 'num2str(nu1)')

dotProd
% save('C:\Users\Nabila Yeazdani\Desktop\trial and error_1');

% toc
find(deltaStore_out>0.3);
