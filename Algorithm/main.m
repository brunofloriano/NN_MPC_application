% This code implements the neural network-based (NN) Model Predictive
% Control (MPC) algorithm for vehicles that enter and exit formation 
% by Bruno R. O. Floriano

clear; close all; clc;
addpath 'functions'

% Receives matrices A,B,L,Pi_estimated,K,
% initial conditions alpha, phi e phi_dot
% and the parameters Delta, delta, tmax, mu e tau

% Uncomment the desired system 
addpath 'Quadrotor'
%addpath 'Car'
%addpath 'GL20_Fig3'
%addpath 'ZH15'


start
save_folder = 'results\';
%%%%%%%%%%%%%%%%%%% CHOOSE SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
horizon = 100;
Kdelta = 24/100; %0.05;


MULTIPLE_SIMULATIONS = 0;
MAX_SIMULATIONS = 20; %if MULTIPLE_SIMULATIONS is on, set the max simulations

% Select the parameter you want to sweep
%PARAMETER_SWEEP = 'horizon = simulation_counter;'; % For horizon analysis
PARAMETER_SWEEP = 'nNuerons = simulation_counter;'; % For cross-validation
%PARAMETER_SWEEP = 'Kdelta = simulation_counter/100;'; % For Kdata analysis

SAVE_DATA = 0;
ONLINE = 1;
MPC = 1;


%%%%%%%%%%%%%%%%%%% INDIVIDUAL SYSTEM DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_state = initial_conditions();
n = size(initial_state{1},1);
m = size(B,2);


%%%%%%%%%%%%%%%%%%%%% COMMUNICATION GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(initial_state,2)-1;
S = length(Pi_estimated);

Psi = transition_matrix(Pi_estimated,delta,Delta);
mc = dtmc(Psi);
mode = simulate(mc,(tmax)/Delta);

Pi_estimated2 = redistribute_weights(Pi_estimated);
Psi2 = transition_matrix(Pi_estimated2,delta,Delta);
mc2 = dtmc(Psi);
mode2 = simulate(mc,(tmax)/Delta);

%%%%%%%%%%%%%%%%%%%% CONTROL SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = zeros(m,n);

%%%%%%%%%%%%%%%%%%%% AUGMENTED SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_state_augmented = [];
for agent_counter = 1+1:N+1
initial_state_augmented = [initial_state_augmented;initial_state{agent_counter}];
end
Aaug = kron(eye(N),A);
Baug = kron(eye(N),B);
desired_state = stack_states(desired_individual_state);

%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if delta_t < Delta
    fprintf('Sampling time has to be greater than Delta \n','s');
    fprintf('Making tdelta = Delta \n','s');
    delta_t = Delta;
end
t = 0:delta_t:tmax-delta_t;
state(:,1) = initial_state_augmented;


if MULTIPLE_SIMULATIONS == 1 % for offline training
    max_sim = MAX_SIMULATIONS;
else
    max_sim = 1;
end

% ONLINE TRAINING LOOP
for simulation_counter = 1:max_sim
    eval(PARAMETER_SWEEP);
    
    text = [t tmax:delta_t:(tmax+horizon*delta_t)];
    clear Jdata net
    Kdata = zeros(1,m*n);
    
    % DATA INITIALIZATION
    state(:,1) = initial_state_augmented;
    individual_state = initial_state;
    state_error(:,1) = state(:,1) - (kron(ones(N,1),individual_state{1}) + desired_state);
    state_error_norm(:,1) = norm(state_error(:,1));
    
    delay = mu/2*randn(length(t),1) + (tau + mu/2);
    theta = mode(round(1*delta_t/Delta));
    Laplacian = L{theta};
    Adjacency = Laplacian2adj(Laplacian);
    J(1) = 0.5*state(:,1)'*kron(Laplacian+Laplacian',eye(n))*state(:,1);
    Jdata = J(1);
    Kline = Kdata;
    Khold = K;
    Jmin = -1e20;
    Jprev = J(1);
    counter = 2;

 
    % TIME LOOP
    for time_counter = 2:length(t)
        
        % ONLINE LEARNING
        if ONLINE == 1
            if time_counter <= 100
                condition = true;
            else
                condition = false; %J(i-1) > 1e4; %mean(abs(J(i-100:i-1))) > 10; %J(i-1) >= -1e30
                vf = 0;
                K(2,:) = zeros(1,n);
            end
            
            if  condition
                net = feedforwardnet(nNeurons);
                net.trainParam.showWindow = 0;
                net = train(net,Kdata',Jdata);
                for k = 1:m*n
                    Kvector2(k,:) = (Kline(k)-Kdelta):Kdelta/2:(Kdelta+Kline(k));
                end
                Ktest = Kvector2(1,:);
                for k = 2:m*n
                    Ktest = combvec(Ktest,Kvector2(k,:));
                end
                y = net(Ktest);
                [value,index] = min(abs(y));
                Kline = Ktest(:,index)';
                K = vector2matrix(Kline,m,n);
                Kdata(counter,:) = Kline;
                Jdata(counter) = J(time_counter-1);

                Khold = vector2matrix(Kdata(counter-1,:),m,n);
                Jmin = J(time_counter-1);
                counter = counter + 1;
            else
                %K = Khold;
            end
        end

        %K=ones(m,n);
        %%%%%%%%%%%%%%%%% DYNAMIC SYSTEM %%%%%%%%%%%%%%%%%%%%%
        if t(time_counter) >= t_change
            theta = mode2(round(time_counter*delta_t/Delta));
        else
            theta = mode(round(time_counter*delta_t/Delta));
            individual_state(N+1) = initial_state(N+1);
        end
        
        Laplacian = L{theta};
        Adjacency = Laplacian2adj(Laplacian);
        
        if ONLINE == 0
            K = Km{theta};
        end
        
        individual_state = group_dynamics(individual_state,K,Adjacency,delta_t);
        
        
        state(:,time_counter) = stack_states(individual_state);
        state_error(:,time_counter) = state(:,time_counter) - (kron(ones(N,1),individual_state{1}) + desired_state);
        state_error_norm(:,time_counter) = norm(state_error(:,time_counter));
        %%%%%%%%%%%%%%%%%%%%%% MPC %%%%%%%%%%%%%%%%%%%%%
        if MPC == 1
            eye_matrix = eye(S);
            virtual_individual_state = individual_state;
            virtual_state(:,1) = state(:,time_counter);
            virtual_state0(:,1) = state(1:n,time_counter);
            vmode = simulate(mc,horizon,'X0',eye_matrix(theta,:));
            %vmode(1) = theta;
            for virtual_time_counter = 1:horizon
                for s = 1:S
                    if virtual_time_counter > 1
                        virtual_Adjacency = Laplacian2adj(L{s});
                        
                        virtual_individual_state = group_dynamics(virtual_individual_state,K,virtual_Adjacency,delta_t);
                        virtual_state(:,virtual_time_counter) = stack_states(virtual_individual_state);
                        virtual_state0(:,virtual_time_counter) = virtual_individual_state{1};
                        
                    end
                    
                    virtual_state_error = virtual_state(:,virtual_time_counter) - (kron(ones(N,1),virtual_state0(:,virtual_time_counter)) + desired_state);
                    Jinst(virtual_time_counter,s) = cost_function(virtual_state_error,kron((L{s}+L{s}'),eye(n)));
                end
                
                if virtual_time_counter == 1
                    Jmean(virtual_time_counter) = Jinst(virtual_time_counter,theta);
                else
                    Jmean(virtual_time_counter) = sum(Jinst(virtual_time_counter,:).*Psi(vmode(virtual_time_counter-1),:));
                end
                
            end
            
            J(time_counter) = sum(Jmean);
            
        else
            J(time_counter) = state(:,time_counter)'*kron(Laplacian,eye(n))*state(:,time_counter);
        end
        
        Jprev = J(time_counter-1);
        
        % Simulation situation display
        current_sim = time_counter/length(t)*100;
        if rem(time_counter,100) == 0 && MULTIPLE_SIMULATIONS == 0
            fprintf('Sim is at %d %% \n',round(current_sim));
        end

    end


    current_sim = simulation_counter;
    if MULTIPLE_SIMULATIONS == 1
        error = state_error_norm;
        squared_error = error.^2;
        
        MSE(:,simulation_counter) =  sum( squared_error' )/length(squared_error) ;
        for st = 1:6
            info = stepinfo(state(st,:),t);
            settling_time(st,simulation_counter) = info.SettlingTime;
        end
        %MSEn(simulation_counter) = norm(MSE(:,simulation_counter));
        
        Jfinal(simulation_counter) = state(:,time_counter)'*state(:,time_counter);
        fprintf('Sim is at %d %% \n',round(100*current_sim/max_sim));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
title('Without offset added')
%plot(t,state([1:2:N*2-1],:)+alpha');
% charc = [':-.'];
% color = ['rbg'];
% vector = [1:2:N*2-1];
% for agent = 1:N
%     plot(t,state(vector(agent),:),[color(agent) charc(agent)]);
%     hold on
% end

%plot(t,state([1:2:N*2-1],:));
plot(t,state([1:n:n*(N-1)+1],:));

xlabel('t')
ylabel('p')
subplot(2,1,2)

%plot(t,state([2:2:N*2],:));
plot(t,state([2:n:n*(N-1)+2],:));

ylabel('v')
xlabel('t')
legend('a1','a2','a3')

figure
semilogy(J)

% figure
% subplot(2,1,1)
% title('With offset added')
% plot(t,state([1:2:N*2-1],:)+alpha');
% %plot(t,state([1:2:N*2-1],:));
% xlabel('t')
% ylabel('\phi_i')
% subplot(2,1,2)
% plot(t,state([2:2:N*2],:));
% ylabel('\phi`_i')
% xlabel('t')

if MULTIPLE_SIMULATIONS == 1
    figure
    semilogy([1:max_sim],MSE)
    xlabel('M')
    ylabel('MSE')
end

% figure
% plot(0:Delta:tmax,mode)
% ylim([0 3])
% xlabel('t')
% ylabel('Markov mode')
% 
% figure
% plot(t,delay)
% xlabel('t')
% ylabel('Delay')

%sim('WW20sim_aug');


%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SAVE_DATA == 1
    comment = ['This simulation was perform with:' newline];
    comment = [comment '']; %Include here what is different for this sim
    clocktime = clock;
    time_now = [num2str(clocktime(1)) '-' num2str(clocktime(2)) '-' num2str(clocktime(3)) '-' num2str(clocktime(4)) '-' num2str(clocktime(5))];
    savefolderlocal = [save_folder '\' save_mainname time_now];
    save_file = [savefolderlocal '\' save_mainname  time_now '.mat'];

    mkdir(savefolderlocal);
    save(save_file);
    saveas(gcf,[savefolderlocal '\' save_mainname  time_now '.eps']);
    saveas(gcf,[savefolderlocal '\' save_mainname  time_now '.fig']);
end
