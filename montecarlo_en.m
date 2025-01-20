clear; clc; close all;

% Define the positions of base stations and users
bs1_pos = [0, 0];      % Position of base station 1 (m)
bs2_pos = [1200, 0];    % Position of base station 2 (m)
radius = 700;          % Coverage radius of the base station (m)
numUsersPerSlot = 500; % Total number of users
ue1_pos = zeros(numUsersPerSlot, 2); % Initialize user positions
numSlots = 1000;       % Total number of time slots
maxdatarate = 2;       % Maximum data rate
lambda = [0.75, 1.5, 2.5];% Message arrival rates (high, medium, low loads)
packet_size_bits = 6.462406251802891e+06;  % Packet size in bits 
slot_duration = 1;         % Duration of each time slot (1 second)

% Parameters for SINR calculation
P_tx1 = 30;  % Transmission power of base station 1
P_tx2 = 30;  % Transmission power of base station 2
P_s = 0; % Default signal power = 0
P_i = 0; % Default interference power = 0
d0 = 1;     % Reference distance
alpha = 3;  % Path loss exponent
N0 = 1e-21;  % Noise power density
B = 5e6;    % Bandwidth
calc_power = @(dist, PL0, exp, d0) 10^((PL0 + 10 * exp * log10(d0 / dist)) / 10); % Power calculation formula

% Define three methods and load conditions
methods = {'Uncoordinated', 'QLBS', 'P-persistent'};
%p_persistent = 0.7; % Transmission probability for P-persistent method

user_success_probability = zeros(3, length(lambda)); % Average success probability for each method
shannon_stastic=zeros(3, length(lambda));
% Initialize user positions
for i = 1:numUsersPerSlot
    r1 = sqrt(rand()) * radius;    % Uniformly distributed distance
    theta1 = rand() * 2 * pi;      % Uniformly distributed angle
    ue1_pos(i, :) = [bs1_pos(1) + r1 * cos(theta1), bs1_pos(2) + r1 * sin(theta1)];
end

% Iterate over three load conditions
for loadIdx = 1:length(lambda)
    arrivalRate = lambda(loadIdx);
    if(loadIdx==1)
        p_persistent=0.9;
    else 
        if(loadIdx==2)
            p_persistent=0.55;
        else 
            if(loadIdx==3)
            p_persistent=0.55;
            end
       end
    end

    % Randomly generate message queues for both base stations
    queueBS1 = poissrnd(arrivalRate, 1, numSlots);
    queueBS2 = poissrnd(arrivalRate, 1, numSlots);

    % Initialize transmission records for base stations
    activeSlots = zeros(1, 3);
    delay1 = zeros(3, 1);
    delay2 = zeros(3, 1);
    for method = 1:3
        successCount = 0; % Total count of successful transmissions
        queue1 = cumsum(queueBS1); % Cumulative queue length for base station 1
        queue2 = cumsum(queueBS2); % Cumulative queue length for base station 2
        transmitBS1 = zeros(1, numSlots);
        transmitBS2 = zeros(1, numSlots);
        shannon_slot = zeros(1,numSlots);
        for t = 1:numSlots
            if method == 1 % Uncoordinated method
            if queue1(t) >= maxdatarate && queue2(t) >= maxdatarate
                transmitBS1(t) = maxdatarate;
                transmitBS2(t) = maxdatarate;
            elseif queue1(t) >= maxdatarate
                transmitBS1(t) = maxdatarate;
                transmitBS2(t) = queue2(t);
            elseif queue2(t) >= maxdatarate
                transmitBS1(t) = queue1(t);
                transmitBS2(t) = maxdatarate;
            else
                transmitBS1(t) = queue1(t);
                transmitBS2(t) = queue2(t);
            end
            elseif method == 2 % Larger method
                if queue1(t) > 0 && queue2(t) > 0
                    if queue1(t) > queue2(t)
                        transmitBS1(t) = min(queue1(t), maxdatarate);
                        transmitBS2(t) = 0;
                    elseif queue1(t) < queue2(t)
                        transmitBS2(t) = min(queue2(t), maxdatarate);
                        transmitBS1(t) = 0;
                    else
                        if rand > 0.5
                            transmitBS1(t) = min(queue1(t), maxdatarate);
                            transmitBS2(t) = 0;
                        else
                            transmitBS2(t) = min(queue2(t), maxdatarate);
                            transmitBS1(t) = 0;
                        end
                    end
                else
                    if queue1(t) > 0 && queue2(t) == 0
                        transmitBS1(t) = min(queue1(t), maxdatarate);
                        transmitBS2(t) = 0;
                    elseif queue1(t) == 0 && queue2(t) > 0
                        transmitBS2(t) = min(queue2(t), maxdatarate);
                        transmitBS1(t) = 0;
                    else
                        transmitBS1(t) = 0;
                        transmitBS2(t) = 0;
                    end
                end
            elseif method == 3 % P-persistent method
                if rand > p_persistent
                    transmitBS1(t) = 0; 
                else 
                    if queue1(t) > maxdatarate
                        transmitBS1(t) = maxdatarate;
                    else
                        transmitBS1(t) = queue1(t);
                    end
                end

                if rand > p_persistent
                    transmitBS2(t) = 0;
                else 
                    if queue2(t) > maxdatarate
                        transmitBS2(t) = maxdatarate;
                    else
                        transmitBS2(t) = queue2(t);
                    end
                end
            end
                queue1(1,:) = queue1(1,:) - transmitBS1(t); % Update queue for base station 1
                queue2(1,:) = queue2(1,:) - transmitBS2(t); % Update queue for base station 2

                if transmitBS1(t) > 0
                    activeSlots(method) = activeSlots(method)+1;
                end
                % Calculate SINR for users
                shannon_c=zeros(1,numUsersPerSlot);
                for u = 1:numUsersPerSlot
                    % Distance between user and base stations
                    dbs1 = norm(ue1_pos(u, :) - bs1_pos);
                    dbs2 = norm(ue1_pos(u, :) - bs2_pos);

                    % Signal and interference power
                    if transmitBS1(t) > 0 && transmitBS2(t) == 0
                        P_s = calc_power(dbs1, P_tx1, alpha, d0);
                        P_i = 0; % No interference
                    elseif transmitBS1(t) > 0 && transmitBS2(t) > 0
                        P_s = calc_power(dbs1, P_tx1, alpha, d0);
                        P_i = calc_power(dbs2, P_tx2, alpha, d0);
                    else
                        break
                    end
                    % Noise power
                    P_n = N0 * B * rand();
                    SINR = 10 * log10(P_s / (P_i + P_n));
                    shannon_c(u)=B*log2(1+10^(SINR/10)); 
                    % Successful transmission condition
                    if SINR > 5
                        successCount = successCount + 1;
                    end     
                end   
                shannon_slot(t)=sum(shannon_c);
                
            end     
        shannon_stastic(method, loadIdx)=sum(shannon_slot)/(numUsersPerSlot*numSlots);
        user_success_probability(method, loadIdx) = successCount / (activeSlots(method) * numUsersPerSlot);
    end
    end
max_bps=max(shannon_stastic(:));

% Display results
% disp('Probability of success of each method under different loads:.');
% disp(user_success_probability);
% 
% figure;
% bar(user_success_probability');
% set(gca, 'XTickLabel', {'Low Load', 'Medium Load', 'High Load'});
% legend(methods, 'Location', 'northwest');
% xlabel('Load Conditions');
% ylabel('Average Success Probability');
% title('Average Success Probability under Different Load Conditions');
% grid on;
% Shannon formula (non-normalized)
figure;
bar(shannon_stastic');
set(gca, 'XTickLabel', {'Low Load', 'Medium Load', 'High Load'});
legend(methods, 'Location', 'northwest');
xlabel('Load Conditions');
ylabel('throughput(bps)');
title('With SINR Threshold=19dB, Shannon formula');
grid on;
% Shannon formula calculates the average throughput of one user per time slot
% Normalized Shannon formula
figure;
bar(shannon_stastic'./max_bps);
set(gca, 'XTickLabel', {'Low Load', 'Medium Load', 'High Load'});
legend(methods, 'Location', 'northwest');
xlabel('Load Conditions');
ylabel('throughput');
title('With SINR Threshold=5dB, The normalized Shannon formula');
grid on;
% The normalized Shannon formula calculates the average throughput of one user per time slot
loadfig = {lambda(1), lambda(2), lambda(3)};
for loadIdx = 1:length(lambda)%simulation
    arrivalRate = lambda(loadIdx);
    queueBS1 = poissrnd(arrivalRate, 1, numSlots);
    queueBS2 = poissrnd(arrivalRate, 1, numSlots);
    transmitBS1 = zeros(1, numSlots);
    transmitBS2 = zeros(1, numSlots);
    transmitstatusBS1=zeros(3, numSlots);
    transmitstatusBS2=zeros(3, numSlots);
    delay1 = zeros(3, 1);
    delay2 = zeros(3, 1);
    for method = 1:3
        queue1 = cumsum(queueBS1);
        queue2 = cumsum(queueBS2);
        for t = 1:numSlots
            if method == 1 % Uncoordinated
                if queue1(t) >= maxdatarate && queue2(t) >= maxdatarate
                    transmitBS1(t) = maxdatarate;
                    transmitBS2(t) = maxdatarate;
                elseif queue1(t) >= maxdatarate
                    transmitBS1(t) = maxdatarate;
                    transmitBS2(t) = queue2(t);  
                elseif queue2(t) >= maxdatarate
                    transmitBS1(t) = queue1(t); 
                    transmitBS2(t) = maxdatarate;
                else
                    transmitBS1(t) = queue1(t);
                    transmitBS2(t) = queue2(t);  
                end
            elseif method == 2 
                if queue1(t) > 0 && queue2(t) > 0
                    if queue1(t) > queue2(t) 
                        transmitBS1(t) = min(queue1(t), maxdatarate);
                        transmitBS2(t)=0;
                    elseif queue1(t) < queue2(t)
                        transmitBS2(t) = min(queue2(t), maxdatarate);
                        transmitBS1(t)=0;
                    else
                        if rand > 0.5
                            transmitBS1(t) = min(queue1(t), maxdatarate);
                            transmitBS2(t)=0;
                        else
                            transmitBS2(t) = min(queue2(t), maxdatarate);
                            transmitBS1(t)=0;
                        end
                    end
                else
                    if queue1(t)>0 && queue2(t)==0
                        transmitBS1(t) = min(queue1(t), maxdatarate);
                        transmitBS2(t)=0;
                    elseif queue1(t)==0 && queue2(t)>0
                        transmitBS2(t) = min(queue2(t), maxdatarate);
                        transmitBS1(t)=0;
                    else
                        transmitBS1(t)=0;
                        transmitBS2(t)=0;
                    end
                end
            elseif method == 3
                if rand > p_persistent
                    transmitBS1(t) = 0; 
                else 
                    if queue1(t) > maxdatarate
                        transmitBS1(t) = maxdatarate;
                    else
                        transmitBS1(t) = queue1(t);
                    end
                end
                if rand > p_persistent
                    transmitBS2(t) = 0;
                else 
                    if queue2(t) > maxdatarate
                        transmitBS2(t) = maxdatarate;
                    else
                        transmitBS2(t) = queue2(t);
                    end
                end
            end 
            queue1(1,:) = queue1(1,:) - transmitBS1(t);
            queue2(1,:) = queue2(1,:) - transmitBS2(t);
            delay1(method,1) = delay1(method,1)+queue1(t);
            delay2(method,1) = delay2(method,1)+queue1(t);
end
        transmitstatusBS1(method,:)=transmitBS1;
        transmitstatusBS2(method,:)=transmitBS2;
    end
    loadcell1{loadIdx} = transmitstatusBS1;
    loadcell2{loadIdx} = transmitstatusBS2;
    delaycell1{loadIdx} = delay1;
    delaycell2{loadIdx} = delay2;
end
% Calculated actual throughput (bits/s)
actual_throughput = zeros(3, 3);
for loadindex = 1:3
    for methodIndex = 1:3
        load1 = loadcell1{loadindex};
        loadrow1 = load1(methodIndex,:);
        load2 = loadcell2{loadindex};
        loadrow2 = load2(methodIndex,:);
         % converted to bits/s
        actual_throughput(loadindex,methodIndex) = user_success_probability(methodIndex, loadindex).*(sum(loadrow1))/numSlots * packet_size_bits / slot_duration;
    end
end
% Normalized separately for each load condition
max_actual = max(actual_throughput(:));


% Plot the actual throughput
figure;
bar(actual_throughput./max_actual);
set(gca, 'XTickLabel', {'Low Load', 'Medium Load', 'High Load'});
legend(methods, 'Location', 'northwest');
xlabel('Load Conditions');
ylabel('Normalized Throughput');
title('With SINR Threshold=5dB, Normalized Actual Throughput');
grid on;