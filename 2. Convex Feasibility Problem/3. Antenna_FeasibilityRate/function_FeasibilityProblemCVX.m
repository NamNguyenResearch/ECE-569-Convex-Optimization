function [feasible,Wsolution] = function_FeasibilityProblemCVX(H,gammavar)
%Solves the feasibility problem with quality-of-service (QoS)

Kr = size(H,1); %Number of users
N = size(H,2);  %Number of transmit antennas 
D = repmat(eye(N),[1 1 Kr]);


%==========================================================================
%Solve the power minimization using CVX
cvx_begin
cvx_quiet(true); 

variable W(N,Kr) complex;   %Variable for N x Kr beamforming matrix
variable POWER              %Scaling parameter for power constraints
minimize POWER              %Minimize the power

subject to

%SINR constraints (Kr constraints)
for k = 1:Kr
    
    %Channels of the signal intended for user i when it reaches user k
    hkD = zeros(Kr,N);
    for i = 1:Kr
        hkD(i,:) = H(k,:)*D(:,:,i);
    end
    
    imag(hkD(k,:)*W(:,k)) == 0; %Useful link 
    
    %SOCP formulation for the SINR constraint of user k
    real(hkD(k,:)*W(:,k)) >= sqrt(gammavar)*norm([1 hkD(k,:)*W(:,[1:k-1 k+1:Kr])  ]);
end

%Power constraints (L constraints) scaled by the variable betavar
norm(W,'fro') <= POWER;
POWER >= 0; %Power constraints must be positive

cvx_end


%==========================================================================
%Analyze result and output variables
%Both power minimization problem and feasibility problem are infeasible
if isempty(strfind(cvx_status,'Solved')) 
    feasible = false;
    Wsolution = [];
else 
    %Both power minimization problem and feasibility problem are feasible
    feasible = true;
    Wsolution = W;
end