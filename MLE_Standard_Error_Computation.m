
function SE_MLE=MLE_Standard_Error_Computation(Theta_Estimate,A_Parameter,b_Parameter,U,Test_Length)
% this function is used to calculate the standard error when using MLE method

% SE_MLE stores the standard error of MLE method
% column vector Theta_Estimate stores the current ability vector estimate of examinee
% matrix A_Parameter stores the discrimination parameters of all items
% column vector b_Parameter stores the b parameters of all items
% column vector U stores the IDs of items which the current examinee has answered
% Test_Length is current test length


Number_of_Dimensions=length(Theta_Estimate);

Fisher_Test_Information_Matrix=zeros(Number_of_Dimensions,Number_of_Dimensions);
SE_MLE=zeros(Number_of_Dimensions,1);
    
if (Test_Length>0)
    
    Item_Answered_ID=U(1:Test_Length,:);
    A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);
    b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);
    
    % compute the item response functions of the Test_Length items evaluated at Theta_Estimate
    IRFs=1./(1+exp(-A_Parameter_Answered*Theta_Estimate).*exp(b_Parameter_Answered));
    
    P_Q=IRFs.*(1-IRFs);
    for i=1:Test_Length
        Fisher_Test_Information_Matrix=Fisher_Test_Information_Matrix+P_Q(i,1)*((A_Parameter_Answered(i,:))'*(A_Parameter_Answered(i,:)));
    end
    
end

Dispersion_Matrix=inv(Fisher_Test_Information_Matrix);

for k=1:Number_of_Dimensions
    SE_MLE(k,1)=sqrt(Dispersion_Matrix(k,k));
end

end
