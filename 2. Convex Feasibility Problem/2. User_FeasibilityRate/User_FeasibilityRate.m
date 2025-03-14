%==========================================================================
close all;
clear all;

numberAntenna = 3;  %Number of Antenna 

K_seq = 95:1:100;
feasibility_ratio = 95:1:100;
K_index = 0;
gamma_dB = -15; 

%==========================================================================
for numberUser=95:1:100 
    
    K_index = K_index + 1;
    number = 0;
    for test=1:20

        H = []; %initialize H matrix

        for i=1:numberUser
            h = 1/sqrt(2*numberUser)*mvnrnd(zeros(numberAntenna,1),eye(numberAntenna),1)'+1i/sqrt(2*numberUser)*mvnrnd(zeros(numberAntenna,1),eye(numberAntenna),1)';
            H = [H h];
        end

        H = H';

        gamma = db2mag(2*gamma_dB);

        [feasible,Wsolution] = function_FeasibilityProblemCVX(H,gamma);
        number = number + feasible;
    end
    feasibility_ratio(K_index) = number/20.0;
end


%==========================================================================
figure
grid on
plot(K_seq,feasibility_ratio,'b-o','LineWidth',1.5);
xlabel('Number of Antenna');
ylabel('Feasibility Ratio');