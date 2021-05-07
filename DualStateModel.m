function DualStateModel(dt,T,alpha,beta,Io)

%Syntax: DualStateModel(dt,T,alpha,beta,Io)
%
%dt    is the simulation time step (in ms) supposing that no more than a 
%      single closure or opening event can occur in this time period(Markov event);
%T     is the total simulation time (in ms);
%alpha is the opening probability in dt (i.e. Po = alpha*dt);
%beta  is the closing probability in dt (i.e. Pc = beta*dt);
%Io    is the current of the open state (in pA).
%
%Example:   DualStateModel(0.1,1000,0.02,0.04,50); to look at a single event
%           DualStateModel(0.1,100000,0.02,0.04,50); to look at the statistics

N = ceil(T/dt);     %Number of steps
state = zeros(1,N); %Array of channel state at all time-steps
curr_state = 0;     %Initial state is set as closed
t = dt;             %State time
T_array = dt*(1:N); %Array of absolute time
Po = alpha*dt;      %Opening probability
Pc = beta*dt;       %Closure probability

%Compute channel state at all time steps
for i =1:N-1
    if curr_state == 0      
        if rand<=Po         
            curr_state = 1; 
            state(i+1) = 1;
            t = dt;
        else
            curr_state = 0; 
            t = t+dt;
        end
    else                    
        if rand <= Pc
            curr_state = 0;
            t = t+dt;
        else
            curr_state = 1;
            state(i+1) = 1;
            t = t+dt;
        end
    end
end

%Plot the temporal evolution of single channel current
figure;
real_current = Io*(state + (randn(1,N))*0.05); %add normally-distributed random noise
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

%Compute open and close state duration
logi = logical(state);           

open_array = diff(logi);  
open_array = find(open_array);  
open_array = diff(open_array);  
open_array(rem(find(open_array),2)==0)=[]; 

clogi = [1 state ~state(end)];
clogi = logical(clogi);

closed_array = diff(clogi); 
closed_array = find(closed_array);  
closed_array = diff(closed_array);  
closed_array(rem(find(closed_array),2)==0)=[]; 

%Close time histogram
figure;
hC = histogram(closed_array*dt); 
t_c = 0: hC.BinEdges(end)/100 : hC.BinEdges(end);
histogram(closed_array*dt, "Normalization", 'pdf');
xlabel('Duration (ms)')
ylabel('pdf for closed state')
title('Statistics of the closed state', 'FontSize',14)
hold on;
plot(t_c, alpha*exp(-alpha*t_c),'r', 'LineWidth',1.5);
legend('Histogram of simulated duration','Theoretical pdf', 'FontSize',12)

%Open time histogram
figure;
hO = histogram(open_array*dt); 
t_o = 0: hO.BinEdges(end)/100 : hO.BinEdges(end);
histogram(open_array*dt, "Normalization", 'pdf');
xlabel('Duration (ms)')
ylabel('pdf for open state')
title('Statistics of the open state', 'FontSize',14)
hold on;
plot(t_o, beta*exp(-beta*t_o),'r','LineWidth',1.5);
legend('Histogram of simulated duration','Theoretical pdf','FontSize',12)

end