
function EAP_Estimate=EAP_Ability_Estimation_Method(A_Parameter,b_Parameter,U,V,Test_Length,Theta_Draw)
% this function is used to estimate the examinees' ability vectors using "Expected A Posterior"(EAP) method

%  column vector EAP_Estimate returns the EAP ability vector estimate
%  matrix A_Parameter stores the discrimination parameters of all items
%  column vector b_Parameter stores the b parameters of all items
%  array U stores the IDs of items which have been answered by the current examinee 
%  array V stores the responses of the current examinee on the items which have been answered by him/her
%  Test_Length is the current test length
%  matrix Theta_Draw stores a random sample drawn from a given distribution


[Number_of_Draws,Number_of_Dimensions]=size(Theta_Draw);

Item_Answered_ID=U(1:Test_Length,:);
V_Answered=V(1:Test_Length,:);
A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);
b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);

IRFs=1./(1+exp(-Theta_Draw*A_Parameter_Answered').*exp(repmat(b_Parameter_Answered',Number_of_Draws,1)));

V=repmat(V_Answered',Number_of_Draws,1);
L=prod((IRFs.^V).*(1-IRFs).^(1-V),2);         % size: Number_of_Draws*1

Numerator=(sum(Theta_Draw.*repmat(L,1,Number_of_Dimensions),1))';
Denominator=sum(L,1);

EAP_Estimate=Numerator/Denominator;


end

