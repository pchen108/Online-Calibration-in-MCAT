
function IRFs_New=Item_Response_Functions_New(Theta_Vector,Item_Parameter)
% this function is used to compute the item response function values of the
% current new item according to multidimensional item reponse theory model, e.g., M2PLM

%  column vector IRFs_New stores the item response function values of the current new item
%  matrix Theta_Vector stores the theta vectors of the examinees who answer the current new item
%  column vector Item_Parameter stores the item parameters of the current new item 


Number_of_Dimensions=length(Theta_Vector(1,:));
Number_of_Examinees=length(Theta_Vector(:,1));

A_Parameter=Item_Parameter(1:Number_of_Dimensions,:);
b_Parameter=Item_Parameter(Number_of_Dimensions+1,:);

Parameter=repmat(A_Parameter',Number_of_Examinees,1);
Exponent_Part=sum(Theta_Vector.*Parameter,2)-b_Parameter;

IRFs_New=exp(Exponent_Part)./(1+exp(Exponent_Part));


end

