clear; clc; close all;
bs1_pos = [0, 0];      
bs2_pos = [1200, 0];
radius = 700;
ue1_pos = [750, 0];
ue2_pos = [450, 0];
numSlots = 1000;
maxdatarate = 2;
lambda = [0.75, 1.5, 2.5];
figure;
hold on; grid on;
plot(bs1_pos(1), bs1_pos(2), 'ro', 'MarkerSize',10,'LineWidth',2);
plot(bs2_pos(1), bs2_pos(2), 'bo', 'MarkerSize',10,'LineWidth',2);

viscircles(bs1_pos, radius, 'Color','r','LineStyle','--');
viscircles(bs2_pos, radius, 'Color','b','LineStyle','--');

plot(ue1_pos(1), ue1_pos(2), 'k*','MarkerSize',3,'LineWidth',2);
plot(ue2_pos(1), ue2_pos(2), 'k^','MarkerSize',3,'LineWidth',2);
title('Network topology diagram');
xlabel('X-coordinate (m)'); ylabel('Y-coordinate (m)');
legend('Base Station 1', 'Base Station 2', 'User 1', 'User 2');
xlim([-1000 2000]);
ylim([-1200 1200]);
hold off;

methods = {'Uncoordinated', 'QLBS', 'P-persistent'};
loadfig = {lambda(1), lambda(2), lambda(3)};
p_persistent = 0.7;
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
%throughput
throughput = zeros(3, 3);
for loadindex = 1:3
    for methodIndex = 1:3
        load1 = loadcell1{loadindex};
        loadrow1 = load1(methodIndex,:);
        load2 = loadcell2{loadindex};
        loadrow2 = load2(methodIndex,:);
        throughput(loadindex,methodIndex) = (sum(loadrow1))/numSlots;
    end
end

figure;
for loadindex = 1:3
    subplot(1, 3, loadindex); 
    bar(throughput(loadindex, :));
    set(gca, 'XTickLabel', methods); 
    xlabel('Methods');
    ylabel('Throughput');
    title('Load index ', loadfig(loadindex));
    ylim([0, max(throughput(:)) * 1.1]); 
    grid on;
end
sgtitle('Throughput Comparison under Different Loads(for BS1 to User1)'); 

delay_results = zeros(3, 3);
load1= zeros(1, 3);
load2= zeros(1, 3);
for loadindex = 1:3
        load1 = delaycell1{loadindex};
        load2 = delaycell2{loadindex};
        delay_results(loadindex,:) = (load1);
        delay_results(loadindex,:) = delay_results(loadindex,:)./(throughput(loadindex,:).*numSlots);
end
figure;
for loadindex = 1:3
    subplot(1, 3, loadindex); 
    bar(delay_results(loadindex,:));
    set(gca, 'XTickLabel', methods); 
    xlabel('Methods');
    ylabel('Average Delay (slots)');
    title('Load index ', loadfig(loadindex));
    %ylim([0, max(loadindex(:)) * 1.1]); 
    grid on;
end
sgtitle('Delay Comparison under Different Loads(for BS1 to User1)');

P_tx1 = 1;  
P_tx2 = 1; 
dbs1_1 = norm(bs1_pos-ue1_pos);    
dbs1_2 = norm(bs1_pos-ue2_pos);    
dbs2_1 = norm(bs2_pos-ue1_pos);  
dbs2_2 = norm(bs2_pos-ue2_pos);  
d0 = 1;     
alpha = 3;  % Path Loss Exponent
N0 = 1e-9;  % noise spectral density
B = 1e3;    % bandwidth
calc_power = @(dist, PL0, exp,d0) 10^((PL0+10*exp*log10(d0/dist))/10);% w
comparison_results = cell(3, 3);  
for loadIndex = 1:3  
    for methodIndex = 1:3     
        data1 = loadcell1{loadIndex}(methodIndex, :);  
        data2 = loadcell2{loadIndex}(methodIndex, :); 
        for i =1 : length(data2)
            if data1(i)>0 && data2(i)>0
                comparison_result(i) = 1;
            else
                comparison_result(i) = 0;
            end
        end
        comparison_results{loadIndex, methodIndex} = comparison_result;
    end
end
sinr_results = cell(3, 3);  
%for user1
for loadIndex = 1:3  
    for methodIndex = 1:3 
         data1 = loadcell1{loadIndex}(methodIndex, :); 
         data1(data1 > 0) = 1;  
         data1(data1 <= 0) = 0;
         for i = 1:length(data1)
             Pt=calc_power(dbs1_1,P_tx1,alpha,d0);
             Pin=calc_power(dbs2_1,P_tx2,alpha,d0);
             sinr(i)=(Pt*data1(i))/((Pin*comparison_results{loadIndex,methodIndex}(i))+N0*B);
         end
         sinr_results{loadIndex, methodIndex} = sinr;
    end
end










