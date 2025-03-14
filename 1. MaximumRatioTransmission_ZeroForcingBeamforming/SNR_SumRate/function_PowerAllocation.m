function powerallocation = function_PowerAllocation(rhos,q,weights)
%Calculates the power allocation 

%Number of base stations
Kt = size(rhos,1); 

%Number of users 
Kr = size(rhos,2); 

%Power allocation
powerallocation=size(Kt,Kr);


%==========================================================================
for j = 1:Kt
    indicesOfNonzero = find(rhos(j,:)>0); %Find users that are served by BS 
    
    %Case 1: Compute waterlevel if all of the users served by BS 
    nuAllActive = (q(j)+sum(1./rhos(j,indicesOfNonzero)))/sum(weights(indicesOfNonzero));
    
    %Case 2: Compute waterlevel if only a subset of the users served by BS
    nuRangeLower = min(1./(rhos(j,indicesOfNonzero)'.*weights(indicesOfNonzero)));
    nuRangeUpper = max(1./(rhos(j,indicesOfNonzero)'.*weights(indicesOfNonzero)));
    nu = fminbnd(@(x) function_AllocDiff(x,q(j),rhos(j,indicesOfNonzero)',weights(indicesOfNonzero)),nuRangeLower,nuRangeUpper);
   
    if function_AllocDiff(nu,q(j),rhos(j,indicesOfNonzero)',weights(indicesOfNonzero)) < function_AllocDiff(nuAllActive,q(j),rhos(j,indicesOfNonzero)',weights(indicesOfNonzero))
        %Compute power allocation with optimal waterlevel (only a subset of users are active)
        powerallocation(j,indicesOfNonzero) = max([weights(indicesOfNonzero)*nu-1./rhos(j,indicesOfNonzero)' zeros(length(indicesOfNonzero),1)],[],2);
    else
        %Compute power allocation with optimal waterlevel (all users are active)
        powerallocation(j,indicesOfNonzero) = max([weights(indicesOfNonzero)*nuAllActive-1./rhos(j,indicesOfNonzero)' zeros(length(indicesOfNonzero),1)],[],2);
    end
    
    %Scale the power allocation 
    powerallocation(j,:) = q(j)*powerallocation(j,:)/sum(powerallocation(j,:));
end


%==========================================================================
function difference = function_AllocDiff(nu,q,rhos,weights)
%Computes absolute difference between the total allocated power and the total
%available powe

difference = abs( sum( max([nu*weights-1./rhos zeros(size(weights))],[],2) ) - q);