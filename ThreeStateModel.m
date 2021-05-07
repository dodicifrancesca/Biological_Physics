function ThreeStateModel(dt,T,alpha,beta,delta,gamma,s,Io)

%Syntax: ThreeStateModel(dt,T,alpha,beta,delta,gamma,s,Io)
%
%dt    is the simulation time step (in ms) supposing that no more than a 
%      single closure or opening event can occur in this time period (Markov event);
%T     is the total simulation time (in ms);
%alpha is the transition rate C->S in dt;
%beta  is the transition rate S->C in dt;
%gamma is the transition rate S->O in dt;
%delta is the transition rate O->S in dt;
%Io    is the current of the fully open state (in pA);
%s     fraction of total current Io in the substate.
%
%Example:   ThreeStateModel(0.1,1000,0.02,0.04,0.03,0.02,0.5,50); to look at a single event
%           ThreeStateModel(0.1,100000,0.02,0.04,0.03,0.02,0.5,50); to look at the
%           statistics

N = ceil(T/dt);         %Number of steps
stateS = zeros(1,N);    %Array of state S (1=channel is in S, 0= channel either C or O)
stateO = zeros(1,N);    %Array of state O (1=channel is in O, 0= channel either C or S)
curr_state = 0;         %Initial state is set as closed (0=C, 1=S, 2=O)
t = dt;                 %State time
T_array = dt*(1:N);     %array of absolute time
Pcs = alpha*dt;         %C->S probability  
Psc = beta*dt;          %S->C probability
Pso = gamma*dt;         %S->O probability
Pos = delta*dt;         %O->S probability
r = rand(1,N);          %Array of random numbers between 0 and 1

%Compute channel state at all time steps
for i =1:N-1
    if curr_state == 0      
        if r(i)<=Pcs        
            curr_state = 1; 
            stateS(i+1) = 1; 
            t = dt;         
        else
            curr_state = 0;  
            t = t+dt;       
        end
    elseif curr_state == 1 
        if r(i)<=Psc       
            curr_state = 0; 
            t = dt;         
        elseif (r(i)>Psc) && (r(i)<=(Psc+Pso))
            curr_state = 2; 
            stateO(i+1) = 1;  
            t = dt;         
        else
            curr_state = 1; 
            stateS(i+1) = 1; 
            t = t+dt;      
        end
    else
        if r(i)<=Pos
            curr_state = 1; 
            stateS(i+1) = 1; 
            t = dt;         
        else
            curr_state = 2; 
            stateO(i+1) = 1; 
            t = t+dt;                    
        end
    
    end
end

state = stateO + s*stateS; %state is 0 when closed, s when S and 1 when open

%Plot the temporal evolution of single channel current
figure;
real_current = Io*(state+(randn(1,N))*0.05); %add normally-distributed random noise
horO = yline(Io,'--r','Open');
horO.LineWidth = 1;
horO = yline(s*Io,'--r','Substate');
horO.LineWidth = 1;
hold on;
horC = yline(0,'--r','Closed');
horC.LineWidth = 1;
hold on;
plot(T_array,real_current);
xlabel('Time (ms)')
ylabel ('Current (pA)')
title ('Single channel current', 'FontSize',14)
axis([-T/10 T+T/3 -Io/2 Io+Io/2]);
pbaspect([2 1 1])
grid on

%Plot statistics
figure;
vertO = xline(Io,'--r','(O)');
vertO.LineWidth = 1;
vertO.LabelOrientation = 'horizontal';
hold on;
vertS = xline(Io*s,'--r','(S)');
vertS.LineWidth = 1;
vertS.LabelOrientation = 'horizontal';
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
%Duration of fully open state
logiO = logical(stateO);   

O_array = diff(logiO); 
O_array = find(O_array); 
O_array = diff(O_array);  
O_array (rem(find(O_array),2)==0)=[];

%Duration of open substate
logiS = logical(stateS);   

S_array = diff(logiS);  
S_array = find(S_array);  
S_array = diff(S_array);  
S_array (rem(find(S_array),2)==0)=[]; 

%Duration of closed state
clogi = [1 state ~state(end)];
clogi = logical(clogi);

C_array = diff(clogi);  
C_array = find(C_array); 
C_array = diff(C_array);  
C_array(rem(find(C_array),2)==0)=[]; 

%Eigenvalue problem
Q = [-alpha beta 0; alpha -beta-gamma delta; 0 gamma -delta];
l = eig(Q);
disp(l);

%Closed time histogram
figure;
hC = histogram(C_array*dt,40);
[f_C, edgC] = histcounts(C_array*dt,40); 
tc = edgC(1:end-1) + diff(edgC)/ 2; 
histogram(C_array*dt,40);
xlabel('Duration (ms)')
ylabel('Absolute frequency')
title('Statistics of the closed state', 'FontSize',14)
hold on;
% Simulated points
tc = tc.';
yc = f_C;
yc = yc.';
plot( tc , yc, 'ro')
hold on;
% Fit with exponential
Fc = @(xc,xdata)xc(1)*exp(l(1)*xdata) + xc(2)*exp(l(1)*xdata) + xc(3)*exp(l(3)*xdata);
x0 = [1 1 1]; %initial conditions
[xc,resnorm,~,~,output] = lsqcurvefit(Fc,x0,tc,yc)

plot(tc,Fc(xc,tc), 'color', 'r','LineWidth',1.5)
legend('Histogram of simulated duration', 'Simulated frequencies' ,'Distribution fit', 'FontSize',12)

%Open substate time histogram
figure;
hS = histogram(S_array*dt,40);
[f_S, edgS] = histcounts(S_array*dt,40); 
ts = edgS(1:end-1) + diff(edgS)/ 2;   
histogram(S_array*dt,40);
xlabel('Duration (ms)')
ylabel('Absolute frequency')
title('Statistics of the substate state', 'FontSize',14)
hold on;
% Simulated points
ts = ts.';
yo = f_S;
yo = yo.';
plot( ts , yo, 'ro')
hold on;
% Fit with exponential
Fs = @(xs,xdata)xs(1)*exp(l(1)*xdata) + xs(2)*exp(l(1)*xdata) + xs(3)*exp(l(3)*xdata);
x0 = [1 1 1]; %initial conditions
[xs,resnorm,~,exitflag,output] = lsqcurvefit(Fs,x0,ts,yo)

plot(ts,Fs(xs,ts), 'color', 'r','LineWidth',1.5)
legend('Histogram of simulated duration', 'Simulated frequencies' ,'Distribution fit', 'FontSize',12)
 
%Open time histogram
figure;
hO = histogram(O_array*dt,40);
[f_O, edgO] = histcounts(O_array*dt,40); 
to = edgO(1:end-1) + diff(edgO)/ 2;   
histogram(O_array*dt,40);
xlabel('Duration (ms)')
ylabel('Absolute frequency')
title('Statistics of the open state', 'FontSize',14)
hold on;
% Simulated points
to = to.';
yo = f_O;
yo = yo.';
plot( to , yo, 'ro')
hold on;
% Fit with exponential
Fo = @(xo,xdata)xo(1)*exp(l(1)*xdata) + xo(2)*exp(l(1)*xdata) + xo(3)*exp(l(3)*xdata);
x0 = [1 1 1]; %initial conditions
[xo,resnorm,~,exitflag,output] = lsqcurvefit(Fo,x0,to,yo) %run solver

plot(to,Fs(xo,to), 'color', 'r','LineWidth',1.5)
legend('Histogram of simulated duration', 'Simulated frequencies' ,'Distribution fit', 'FontSize',12)

end