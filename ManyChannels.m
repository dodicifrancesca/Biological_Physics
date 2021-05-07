function ManyChannels(dt,T,alpha,beta,gamma,delta,C,Io)

%Syntax: ManyChannels(dt,T,alpha,beta,gamma,delta,C,Io)
%
%dt    is the simulation time step (in ms) supposing that no more than a 
%      single closure or opening event can occur in this time period (Markov event);
%T     is the total simulation time (in ms);
%alpha is the transition rate C->O1 in dt;
%beta  is the transition rate O1->C in dt;
%gamma is the transition rate O1->O2 in dt;
%delta is the transition rate O2->O1 in dt;
%C     is the number of channel simulated
%Io    is the current of the open state (in pA).
%
%Example:   ManyChannels(0.1,1000,0.02,0.04,0.002,0.004,5,50); to look at 5 channels
          
N = ceil(T/dt);         %Number of steps
states = zeros(C,N);    %Array of states (1=channel is in O1, 0= channel either C or O2)
curr_state = 0;         %state at current step (0=C, 1=O1, 2=O2)
t = dt;                 %State time
T_array = dt*(1:N);     %Array of absolute time
Pc1 = alpha*dt;         %C -> O1 probability  
P1c = beta*dt;          %O1-> C  probability
P12 = gamma*dt;         %O1-> O2 probability
P21 = delta*dt;         %O2-> O1 probability

%Compute channel state at all time steps
for rep = 1:C-1
    r = rand(1,N);          %array of random numbers between 0 and 1
    for i =1:N-1
        if curr_state == 0      
            if r(i)<=Pc1        
                curr_state = 1; 
                states(rep,i+1) = 1; 
                t = dt;         
            else
                curr_state = 0; 
                t = t+dt;       
            end
        elseif curr_state == 1 
            if r(i)<=P1c       
                curr_state = 0; 
                t = dt;         
            elseif (r(i)>P1c) && (r(i)<=(P1c+P12))
                curr_state = 2; 
                states(rep,i+1) = 1;  
                t = dt;         
            else
                curr_state = 1; 
                states(rep,i+1) = 1;  
                t = t+dt;       
            end
        else
            if r(i)<=P21
                curr_state = 1; 
                states(rep,i+1) = 1; 
                t = dt;         
            else
                curr_state = 2; 
                states(rep,i+1) = 1; 
                t = t+dt;                    
            end    
        end
    end
end

current = sum(states);
real_current = Io*current + Io*(randn(1,N))*0.05; %add normally-distributed random noise
real_current = real_current(2:end); %remove 0 at beginning 
T_array = T_array(2:end);
Iavg = mean(real_current);

%Plot the temporal evolution of the total current
figure;
plot(T_array,real_current);
horO = yline(Iavg,'--r',['Average'; 'current']);
horO.LineWidth = 1;
hold on;
xlabel('Time (ms)')
ylabel ('Current (pA)')
title (['Current from ', num2str(C) ,' channels'], 'FontSize',14)
xlim([-T/10 T+T/5]);
pbaspect([2 1 1])
grid on


end