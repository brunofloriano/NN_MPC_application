% This code implements the neural network-based (NN) Model Predictive
% Control (MPC) algorithm for vehicles that enter and exit formation 
% by Bruno R. O. Floriano

clear; close all; clc;
addpath functions
addpath functions/Bens_functions/
t_change = 1e100;

%file_name = 'dataBalloon 2023-4-24-22-45';
file_name = 'dataBalloon 2023-4-27-3-15';
file = ['results/' file_name '/' file_name '.mat'];

load(file);

% Receives matrices A,B,L,Pi_estimated,K,
% initial conditions alpha, phi e phi_dot
% and the parameters Delta, delta, tmax, mu e tau

% Uncomment the desired system 
%addpath 'Quadrotor'
addpath 'Balloon'
%addpath 'Car'
%addpath 'GL20_Fig3'
%addpath 'ZH15'

start_continue

t = 0:delta_t:tmax;%-delta_t;
Time.t = t;

%% COMMUNICATION GRAPH 

%N = size(initial_state,2)-1;
S = length(Pi_estimated);

Psi = transition_matrix(Pi_estimated,delta,Delta);
mc = dtmc(Psi);
mode = simulate(mc,(tmax)/Delta);

Pi_estimated2 = redistribute_weights(Pi_estimated);
Psi2 = transition_matrix(Pi_estimated2,delta,Delta);
mc2 = dtmc(Psi);
mode2 = simulate(mc,(tmax)/Delta);


% ONLINE TRAINING LOOP
for simulation_counter = 1:sim_step:max_sim
 
    % TIME LOOP
    for time_counter = 502:length(t)
        Time.time_counter = time_counter;
        % ONLINE LEARNING
        if ONLINE == 1
            if time_counter <= 100
                condition = true;
            else
                condition = J(time_counter-1) > 10; %true;% false; %J(i-1) > 1e4; %mean(abs(J(i-100:i-1))) > 10; %J(i-1) >= -1e30
                vf = 0;
                %K = K*0;
                %K(2,:) = zeros(1,n);
            end
            
            if  condition
                net = feedforwardnet(nNeurons);
                net.trainParam.showWindow = 0;
                net = train(net,Kdata',Jdata);
                for k = 1:m*nr
                    %Kvector2(k,:) = (Kline(k)-Kdelta):Kdelta/2:(Kdelta+Kline(k));
                    Kvector2(k,:) = (Kline(k)-Kdelta):Kdelta:(Kdelta+Kline(k));
                end
                Ktest = Kvector2(1,:);
                for k = 2:m*nr
                    Ktest = combvec(Ktest,Kvector2(k,:));
                end
                y = net(Ktest);
                [value,index] = min(abs(y));
                Kline = Ktest(:,index)';
                K = vector2matrix(Kline,m,nr);
                Kdata(counter,:) = Kline;
                Jdata(counter) = J(time_counter-1);

                Khold = vector2matrix(Kdata(counter-1,:),m,nr);
                Jmin = J(time_counter-1);
                counter = counter + 1;
            else
                %K = Khold;
            end
        end

        %K=ones(m,n);
        %%%%%%%%%%%%%%%%% DYNAMIC SYSTEM %%%%%%%%%%%%%%%%%%%%%
        if t(time_counter) >= t_change
            theta = mode2(round((time_counter-1)*delta_t/Delta));
        else
            theta = mode(round((time_counter-1)*delta_t/Delta));
            individual_state(N+1) = initial_state(N+1);
        end
        
       
        if strcmp(save_mainname,'dataBalloon ')
            Laplacian = laplacian(Ggraph{time_counter-1});
            Laplacian_time{time_counter} = Laplacian;
            L_cumulative = accumulate_laplacian(Laplacian,L_cumulative);
            L_cumulative_time{time_counter} = L_cumulative;
            K(1) = 0;
            %K = K*0;
            if t(time_counter) > t_change
                L_cumulative = zeros(N);
                encompassed_area = zeros(map_steps);
                t_change = 1e10;
            end
        else
            Laplacian = L{theta};
        end
        
        Adjacency = Laplacian2adj(Laplacian);
        
        if ONLINE == 0
            K = Km{theta};
        end
        
        individual_state = group_dynamics(individual_state,K,Adjacency,delta_t);
        %balloon_kinematics
        [individual_coord, individual_angular_position, individual_state, Ggraph] = group_kinematics(individual_angular_position, individual_state, individual_coord, time_counter, delta_t, comm_range, max_root_connections,Ggraph);
        
        if time_counter == 502
            Ggraph{502} = routing_protocol(position,0,max_root_connections);
            L_cumulative = L_cumulative*0;
        end

        state(:,time_counter) = stack_states(individual_state);
        state_error(:,time_counter) = state(:,time_counter) - (kron(ones(N,1),individual_state{1}) + desired_state);
        state_error_norm(:,time_counter) = norm(state_error(:,time_counter));

        if strcmp(save_mainname,'dataBalloon ')
            %energy_time = time_energy(alpha_time_constant, t(time_counter), time_horizon);
            Area = area_coverage(G,individual_coord, Area, Time, Agents);
            %encompassed_area = encompassed_area*energy_time;
            encompassed_area_time{time_counter} = Area.encompassed_area;
            encompassed_area = Area.encompassed_area;
        end
        %%%%%%%%%%%%%%%%%%%%%% MPC %%%%%%%%%%%%%%%%%%%%%
        if MPC == 1
            eye_matrix = eye(S);
            virtual_individual_state = individual_state;
            virtual_state(:,1) = state(:,time_counter);
            virtual_state0(:,1) = state(1:n,time_counter);
            vmode = simulate(mc,horizon,'X0',eye_matrix(theta,:));
            if strcmp(save_mainname,'dataBalloon ')
                virtual_area_covered = Area.area_covered;
                virtual_encompassed_area = encompassed_area;
                virtual_individual_coord = individual_coord;
                virtual_individual_angular_position = individual_angular_position;
                Virtual_Area = Area;
                Virtual_Time = Time;
                virtual_Ggraph = Ggraph;
            end
            %vmode(1) = theta;
            for virtual_time_counter = 1:horizon
                virtual_time = t(time_counter) + (virtual_time_counter - 1)*delta_t;
                Virtual_Time.time_counter = virtual_time_counter;
                for s = 1:S
                    if virtual_time_counter > 1
                        virtual_Adjacency = Laplacian2adj(L{s});
                        
                        virtual_individual_state = group_dynamics(virtual_individual_state,K,virtual_Adjacency,delta_t);
                        if strcmp(save_mainname,'dataBalloon ')
                            %virtual_ballon_kinematics;
                            [virtual_individual_coord, virtual_individual_angular_position, virtual_individual_state, virtual_Ggraph] = group_kinematics(virtual_individual_angular_position, virtual_individual_state, virtual_individual_coord, virtual_time_counter, delta_t, comm_range, max_root_connections,virtual_Ggraph);
                            %[virtual_area_covered, virtual_encompassed_area] = area_coverage(G,virtual_individual_coord,max_radius,map_steps, virtual_encompassed_area, area_radius,virtual_time_counter,total_circle_area, alpha_time_constant, virtual_time, time_horizon,delta_t,area_covariance_matrix,max_area_pdf,std_area);
                            Virtual_Area = area_coverage(G,virtual_individual_coord, Virtual_Area, Virtual_Time, Agents);
                            virtual_encompassed_area = Virtual_Area.encompassed_area;
                            virtual_area_covered = Virtual_Area.area_covered;
                            %virtual_encompassed_area = virtual_encompassed_area*energy_time;
                        end
                        virtual_state(:,virtual_time_counter) = stack_states(virtual_individual_state);
                        virtual_state0(:,virtual_time_counter) = virtual_individual_state{1};
                        
                    end
                    
                    virtual_state_error = virtual_state(:,virtual_time_counter) - (kron(ones(N,1),virtual_state0(:,virtual_time_counter)) + desired_state);
                    Jinst(virtual_time_counter,s) = cost_function(virtual_state_error,kron((L{s}+L{s}'),eye(n)));
                    cost_function_balloon;
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
        if rem(time_counter,1) == 0 && MULTIPLE_SIMULATIONS == 0
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

%% PLOT
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
    vector = 1:sim_step:max_sim;
    semilogy(vector,MSE(vector))
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


%% SAVE
if SAVE_DATA == 1 || GRAPH_MOVIE == 1
    comment = ['This simulation was perform with:' newline];
    comment = [comment '']; %Include here what is different for this sim
    clocktime = clock;
    time_now = [num2str(clocktime(1)) '-' num2str(clocktime(2)) '-' num2str(clocktime(3)) '-' num2str(clocktime(4)) '-' num2str(clocktime(5))];
    savefolderlocal = [save_folder '/' save_mainname time_now];
    save_file = [savefolderlocal '/' save_mainname  time_now '.mat'];

    mkdir(savefolderlocal);
    save(save_file);
    saveas(gcf,[savefolderlocal '/' save_mainname  time_now '.eps']);
    saveas(gcf,[savefolderlocal '/' save_mainname  time_now '.fig']);
end

%% Generate movie plot
if GRAPH_MOVIE == 1
    plot_movie(savefolderlocal,VIDEO,GIF,t,individual_coord,Ggraph,GRAPH_AREA,Nr,encompassed_area_time)
end