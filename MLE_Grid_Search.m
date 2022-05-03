
function MLE_Estimate=MLE_Grid_Search(Theta_All_Possible,A_Parameter,b_Parameter,U,V,Test_Length)
% this function is used to find the MLE estimates by using grid search method

% column vector MLE_Estimate stores the MLE estimate of examinee
% matrix Theta_All_Possible stores all possible ability vectors
% matrix A_Parameter stores the discrimination parameters of all items
% column vector b_Parameter stores the b parameters of all items
% column vector U stores the IDs of items which the current examinee has answered
% column vector V records the responses to the items which the current examinee has answered
% Test_Length is current test length


Number_of_Ability_Vectors=length(Theta_All_Possible(1,:));

Item_Answered_ID=U(1:Test_Length,:);
V_Answered=V(1:Test_Length,:);
A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);
b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);

IRFs=1./(1+exp(-Theta_All_Possible'*A_Parameter_Answered').*exp(repmat(b_Parameter_Answered',Number_of_Ability_Vectors,1)));

V=repmat(V_Answered',Number_of_Ability_Vectors,1);
L=prod((IRFs.^V).*(1-IRFs).^(1-V),2);

[~,Maximum_Index]=max(L);
MLE_Estimate=Theta_All_Possible(:,Maximum_Index);


end
