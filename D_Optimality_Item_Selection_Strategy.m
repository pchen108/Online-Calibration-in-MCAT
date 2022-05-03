
function Item_ID=D_Optimality_Item_Selection_Strategy(Item_Flag,Theta_Estimate,A_Parameter,b_Parameter,Sigma,U,Test_Length)
% this function is employed to select the next item by using the D-Optimality item selection rule

% Item_ID returns the ID of the selected item
% column vector Item_Flag indicates which items have been answered and which have not
% column vector Theta_Estimate stores the estimated theta vector of examinee
% matrix A_Parameter stores all discrimination parameters of all items
% column vector b_Parameter stores all b parameters of all items
% matrix Sigma is the covariance matrix of the prior distribution of theta
% column vector U stores the IDs of items which the current examinee has answered
% Test_Length is current test length


% Step1: calculate the Fisher Test Information Matrix after answering k-1 items
Fisher_Test_Information_Matrix=Fisher_Test_Information_Computation(Theta_Estimate,A_Parameter,b_Parameter,U,Test_Length);

% Step 2: implement the D-Optimality item selection algorithm
Index=find(Item_Flag==0);                             % ID of items which are available for selection
Number_of_Remaining_Items=length(Index);

A_Parameter_Remaining=A_Parameter(Index,:);
b_Parameter_Remaining=b_Parameter(Index,:);

Determinant_Values=zeros(Number_of_Remaining_Items,1);

IRFs_Remaining=1./(1+exp(-A_Parameter_Remaining*Theta_Estimate).*exp(b_Parameter_Remaining));
P_Q_Remaining=IRFs_Remaining.*(1-IRFs_Remaining);

for i=1:Number_of_Remaining_Items
    Fisher_Item_Information_Matrix=P_Q_Remaining(i,1)*((A_Parameter_Remaining(i,:))'*(A_Parameter_Remaining(i,:)));
    Determinant_Values(i,1)=det(Fisher_Test_Information_Matrix+Fisher_Item_Information_Matrix+inv(Sigma));
end

[~,Maximum_DET_Index]=max(Determinant_Values);
Item_ID=Index(Maximum_DET_Index,:);

end


function Fisher_Test_Information_Matrix=Fisher_Test_Information_Computation(Theta_Estimate,A_Parameter,b_Parameter,U,Test_Length)
% this function is used to compute the Fisher test information matrix based on the answered k-1 items

% matrix Fisher_Test_Information_Matrix returns the final result
% column vector Theta_Estimate stores the estimated theta vector of examinee
% matrix A_Parameter stores all discrimination parameters of all items
% column vector b_Parameter stores all B parameters of all items
% column vector U stores the IDs of items which the current examinee has answered
% Test_Length is current test length


Number_of_Dimensions=length(Theta_Estimate);

Fisher_Test_Information_Matrix=zeros(Number_of_Dimensions,Number_of_Dimensions);
    
if (Test_Length>0)
    
    Item_Answered_ID=U(1:Test_Length,:);
    A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);
    b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);
    
    % compute the item response functions of the Test_Length items at Theta_Estimate
    IRFs=1./(1+exp(-A_Parameter_Answered*Theta_Estimate).*exp(b_Parameter_Answered));
    
    P_Q=IRFs.*(1-IRFs);
    for i=1:Test_Length
        Fisher_Test_Information_Matrix=Fisher_Test_Information_Matrix+P_Q(i,1)*((A_Parameter_Answered(i,:))'*(A_Parameter_Answered(i,:)));
    end
    
end

end

