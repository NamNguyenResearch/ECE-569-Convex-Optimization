%==========================================================================
close all;
clear all;

numberUser = 50;    %Number of User
numberAntenna = 3;  %Number of Antenna 

gamma_dB_Seq = -13:0.5:-11;
feasibility_ratio = -13:0.5:-11;
gamma_index = 0;


%==========================================================================
for gamma_dB=-13:0.5:-11
    gamma_index = gamma_index + 1;
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
    feasibility_ratio(gamma_index) = number/20.0;
end

%==========================================================================
figure
grid on
plot(gamma_dB_Seq,feasibility_ratio,'b-o','LineWidth',1.5);
xlabel('Average SNR [dB]');
ylabel('Feasibility Ratio');