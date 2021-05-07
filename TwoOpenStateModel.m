function TwoOpenStateModel(dt,T,alpha,beta,gamma,delta,Io)

%Syntax: TwoOpenStateModel(dt,T,alpha,beta,gamma,delta,Io)
%
%dt    is the simulation time step (in ms) supposing that no more than a 
%      single closure or opening event can occur in this time period (Markov event);
%T     is the total simulation time (in ms);
%alpha is the transition rate C->O1 in dt;
%beta  is the transition rate O1->C in dt;
%gamma is the transition rate O1->O2 in dt;
%delta is the transition rate O2->O1 in dt;
%Io    is the current of the open state (in pA).
%
%Example:   TwoOpenStateModel(0.1,1000,0.02,0.02,0.004,0.004,50); to look at a single event
%           TwoOpenStateModel(0.1,100000,0.02,0.02,0.004,0.004,50); to look at the
%           statistics

N = ceil(T/dt);         %Number of steps
stateO1 = zeros(1,N);   %Array of state O1 (1=channel is in O1, 0= channel either C or O2)
stateO2 = zeros(1,N);   %Array of state O2 (1=channel is in O2, 0= channel either C or O1)
curr_state = 0;         %state at current step (0=C, 1=O1, 2=O2)
t = dt;                 %State time
T_array = dt*(1:N);     %array of absolute time
Pc1 = alpha*dt;         %C -> O1 probability  
P1c = beta*dt;          %O1-> C  probability
P12 = gamma*dt;         %O1-> O2 probability
P21 = delta*dt;         %O2-> O1 probability
r = rand(1,N);          %array of random numbers between 0 and 1

%Compute channel state throughout simulation time
for i =1:N-1
    if curr_state == 0      %if state at current iteration is closed
        if r(i)<=Pc1        
            curr_state = 1; 
            stateO1(i+1) = 1; 
            t = dt;         
        else
            curr_state = 0; 
            t = t+dt;       
        end
    elseif curr_state == 1  %if state at current iteration is O1
        if r(i)<=P1c       
            curr_state = 0; 
            t = dt;         
        elseif (r(i)>P1c) && (r(i)<=(P1c+P12)) 
            curr_state = 2; 
            stateO2(i+1) = 1;
            t = dt;         
        else
            curr_state = 1; 
            stateO1(i+1)= 1;
            t = t+dt;      
        end
    else                    %if state at current iteration is O2
        if r(i)<=P21
            curr_state = 1; 
            stateO1(i+1)= 1;
            t = dt;         
        else
            curr_state = 2; 
            stateO2(i+1)= 1;
            t = t+dt;                    
        end
    
    end
end

stateO = stateO1 + stateO2; %stateO is 0 when closed and 1 when open (O1 or O2)


%Plot the temporal evolution of single channel current
figure;
real_current = Io*(stateO + (randn(1,N))*0.05); %randn to add normally-distributed random noise
horO = yline(Io,'--r','Open');
horO.LineWidth = 1;
hold on;
horC = yline(0,'--r','Closed');
horC.LineWidth = 1;
hold on;
plot(T_array,real_current);
xlabel('Time (ms)')
ylabel ('Current (pA)')
title ('Single channel current', 'FontSize',14)
axis([-T/10 T+T/5 -Io/2 Io+Io/2]);
pbaspect([2 1 1])
grid on

%Plot statistics
figure;
vertO = xline(Io,'--r','(O)');
vertO.LineWidth = 1;
vertO.LabelOrientation = 'horizontal';
hold on;
vertC = xline(0,'--r','(C)');
vertC.LineWidth = 1;
vertC.LabelOrientation = 'horizontal';
hold on;
histogram(real_current, 100);
xlabel('Current (pA)')
ylabel ('Absolute frequency')
title ('Histogram of the channel current', 'FontSize',14)


%Compute duration of the three states
%OPEN TOT
logi = logical(stateO);           

open_array = diff(logi);  
open_array = find(open_array);  
open_array = diff(open_array);  
open_array(rem(find(open_array),2)==0)=[]; 

%CLOSED
clogi = [1 stateO ~stateO(end)]; %to count also first closure event
clogi = logical(clogi);

closed_array = diff(clogi);  %array with 0 btwn steps where state doesn't change, 1 when 0->1, -1 when 1->0
closed_array = find(closed_array);  %array with positions of nonzero elements (where state changes)
closed_array = diff(closed_array);  
closed_array(rem(find(closed_array),2)==0)=[]; %at odd positions we have durations of closed states

%OPEN
state_durat = stateO1 + 3*stateO2;  %now C = 0, O1 = 1, O2 = 3
state_durat = diff(state_durat);    %array where 0: no transition,   
                                    %            1: trans C -> O1,  -1: trans O1 -> C,   
                                    %            2: trans O1 -> O2, -2: trans O2 -> O1
%OPEN1
O1_durat = find(state_durat==1|state_durat==-1);  %an open state starting with O1 goes from 1 to -1
O1_durat = diff(O1_durat);  
O1_durat (rem(find(O1_durat),2)==0)=[]; %at odd positions we have durations of open states

%OPEN2
O2_durat = [];
for i = 1:N-1
    if state_durat(i)==2
        for j = i:N-1
            if state_durat(j)== -1
                a = j-i;
                O2_durat =[O2_durat a];
                i = j+1;
                break
            end
        end
    end
end

%Define coefficients for solution of differential equation
b = beta + delta + gamma;
D = sqrt(b^2 - 4*beta*delta);
lambda1 = (-b+D)/2;
lambda2 = -(b+D)/2;

a1 = -(lambda2 + beta + gamma)/(lambda1 - lambda2);
a2 = (lambda1 + beta + gamma)/(lambda1 - lambda2);
b1 = gamma/(lambda1 - lambda2);
b2 = -gamma/(lambda1 - lambda2);

%Closed time histogram
figure;
[f_C, edgC] = histcounts(closed_array*dt,50); 
t_c = edgC(1:end-1) + diff(edgC)/ 2; 
histogram(closed_array*dt,50,'Normalization','pdf');
xlabel('Duration (ms)')
ylabel('Absolute frequency')
ylabel('pdf for closed state')
title('Statistics of the closed state', 'FontSize',14)
hold on;
plot(t_c, alpha*exp(-alpha*t_c),'r','LineWidth',1.5);
legend('Histogram of simulated duration','Theoretical pdf','FontSize',12)

%O1 time histogram
figure;
[f_O1, edgO1] = histcounts(O1_durat*dt,80,'Normalization','pdf'); 
t_o1 = edgO1(1:end-1) + diff(edgO1)/ 2;     
histogram(O1_durat*dt, 80, 'Normalization','pdf');
xlabel('Duration (ms)')
ylabel('pdf for state O1')
title('Statistics of the O1 state', 'FontSize',14)
hold on;
PO1 = a1*exp(lambda1*t_o1) + a2*exp(lambda2*t_o1);
Area = trapz(t_o1,PO1);
PO1_normalized = PO1./Area;
plot(t_o1, PO1_normalized,'r','LineWidth',1.5);
legend('Histogram of simulated duration','Theoretical pdf', 'FontSize',12)

%O2 time histogram
figure;
[f_O2, edgO2] = histcounts(O2_durat*dt,25,'Normalization','pdf');
t_o2 = edgO2(1:end-1) + diff(edgO2)/ 2;  
histogram(O2_durat*dt,25,'Normalization','pdf');
xlabel('Duration (ms)')
ylabel('pdf for state O2')
title('Statistics of the O2 state', 'FontSize',14)
hold on;
t = linspace(0,edgO2(end-1),100);
PO2 = b1*exp(lambda1*t) + b2*exp(lambda2*t);
Area = trapz(t,PO2);
PO2_normalized = PO2./Area;
plot(t, PO2_normalized,'r','LineWidth',1.5);
legend('Histogram of simulated duration','Theoretical pdf', 'FontSize', 12)

%Total open time histogram
figure;
[f_O, edgO] = histcounts(open_array*dt,50); 
t_o = edgO(1:end-1) + diff(edgO)/ 2; 
histogram(open_array*dt,50,'Normalization','pdf');
xlabel('Duration (ms)')
ylabel('pdf for open state')
title('Statistics of the open state', 'FontSize',14)
hold on;
PO = -(lambda1*(a1*exp(lambda1*t_o) + b1*exp(lambda1*t_o)) + lambda2*(a2*exp(lambda2*t_o) + b2*exp(lambda2*t_o)));
plot(t_o, PO,'r','LineWidth',1.5);
legend('Histogram of simulated duration','Theoretical pdf', 'FontSize',12)

end