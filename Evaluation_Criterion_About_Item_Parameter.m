
function [MSE,Bias,Correlation_Coefficient]=Evaluation_Criterion_About_Item_Parameter(Item_Parameter_Estimate,Item_Parameter_True)


Number_of_Replications=length(Item_Parameter_Estimate);
[Number_of_New_Items, Number_of_Dimensions]=size(Item_Parameter_Estimate{1,1});

MSE=zeros(1,Number_of_Dimensions);
Bias=zeros(1,Number_of_Dimensions);
Correlation_Coefficient=zeros(1,Number_of_Dimensions);

for rep=1:Number_of_Replications
    
    MSE=MSE+(sum((Item_Parameter_Estimate{rep,1}-Item_Parameter_True{rep,1}).^2,1))/Number_of_New_Items;
    Bias=Bias+(sum(Item_Parameter_Estimate{rep,1}-Item_Parameter_True{rep,1},1))/Number_of_New_Items;

    Correlation_Coefficient_Temp=zeros(1,Number_of_Dimensions);
    for k=1:(Number_of_Dimensions)
        Correlation_Coefficient_Matrix=corrcoef(Item_Parameter_Estimate{rep,1}(:,k),Item_Parameter_True{rep,1}(:,k));
        Correlation_Coefficient_Temp(1,k)=Correlation_Coefficient_Matrix(1,2);
    end
    Correlation_Coefficient=Correlation_Coefficient+Correlation_Coefficient_Temp;
    
end

MSE=MSE/Number_of_Replications;
Bias=Bias/Number_of_Replications;
Correlation_Coefficient=Correlation_Coefficient/Number_of_Replications;

end
