%==========================================================================
close all;
clear all;

numberUser = 100; %Number of Antenna

N_seq = 7:1:12;
feasibility_ratio = 7:1:12;
N_index = 0;
gamma_dB = -10; 

%==========================================================================
for numberAntenna=7:1:12 
   
    N_index = N_index + 1;
    number = 0;
    for test=1:20
        
        %Initialize H matrix
        H = []; 

        for i=1:numberUser
            h = 1/sqrt(2*numberUser)*mvnrnd(zeros(numberAntenna,1),eye(numberAntenna),1)'+1i/sqrt(2*numberUser)*mvnrnd(zeros(numberAntenna,1),eye(numberAntenna),1)';
            H = [H h];
        end

        H = H';

        gamma = db2mag(2*gamma_dB);

        [feasible,Wsolution] = function_FeasibilityProblemCVX(H,gamma);
        number = number + feasible;
    end
    feasibility_ratio(N_index) = number/20.0;
end


grid on

%==========================================================================
figure
grid on
plot(N_seq,feasibility_ratio,'b-o','LineWidth',1.5);
xlabel('Number of Antenna');
ylabel('Feasibility Ratio');