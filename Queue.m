%MMk

%%
total_packets = 1000;
servers = 3;
monte_carlo_simulations = 1;

meanlen = [];
min_mu = 1;
max_mu = 30;
step_mu = 0.5;
min_lambda=1;
max_lambda=30;
step_lambda = 0.5;
[lambda,mu] = meshgrid(min_lambda:step_lambda:max_lambda,min_mu:step_mu:max_mu);




i = 0;
j = 0;
for x = min_lambda:step_lambda:max_lambda
    i = i + 1;
    j = 0 ;
    for y = min_mu:step_mu:max_mu
        j = j + 1;
        meanlength_montecarlo = 0;
        for m = 1:monte_carlo_simulations
            [meanlength, len, leng] = MMK(x,y,servers,total_packets);
            %figure
            %stairs(len(:,1),leng,'-*')
            %hold on
            %stairs(len(:,1),len(:,2),'^')
            meanlength_montecarlo=meanlength_montecarlo+meanlength;
        end
        meanlength = meanlength_montecarlo/monte_carlo_simulations;
        [nq_expected,p0, pQ, rho] = expected_mean_queue(1/x, 1/y, servers);
        %[meanlength,len, leng] = MMK(x,y,servers,total_packets);
        meanlen(i,j) = meanlength;

    end
end
    
figure
hold on
grid on
surf(lambda,mu,meanlen');
xlabel('lambda');
ylabel('mu');
zlabel('mean length of the queue in teory');
hold off
figure
surf(lambda,mu,double([1./lambda>servers*1./mu]));
    
% figure
% hold on
% grid on
% surf(lambda,mu,meanlen_teory');
% xlabel('lambda');
% ylabel('mu');
% zlabel('mean length of the queue in teory');


%%
clear all
clc

total_packets = 1000;
servers = 5;
monte_carlo_simulations = 10;

meanlen = [];
meanlen_theory = [];
meanlen_montecarlo = [];

lambda = [10 20 30 40 50 60];
mu = [1:1:60];


i = 0;
j = 0;
for x = lambda
    i = i + 1;
    j = 0 ;
    for y = mu
        j = j + 1;
        meanlength_montecarlo = 0;
        for m = 1:monte_carlo_simulations
            [meanlength, len, leng] = MMK(x,y,servers,total_packets);
            %figure
            %stairs(len(:,1),leng,'-*')
            %hold on
            %stairs(len(:,1),len(:,2),'^')
            meanlength_montecarlo=meanlength_montecarlo+meanlength;
        end
        meanlength = meanlength_montecarlo/monte_carlo_simulations;
        [nq_expected,p0, pQ, rho] = expected_mean_queue(1/x, 1/y, servers);
        %[meanlength,len, leng] = MMK(x,y,servers,total_packets);
        meanlen(i,j) = meanlength;
        if rho <1-0.01
            meanlen_theory(i,j) = nq_expected;
        else
            meanlen_theory(i,j) = 0;
        end

    end
end


figure
semilogy(meanlen_theory','--')
hold on
semilogy(meanlen')



