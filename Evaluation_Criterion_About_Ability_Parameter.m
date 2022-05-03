
function [MSE,Bias,Correlation_Coefficient,Mean_ED]=Evaluation_Criterion_About_Ability_Parameter(Theta_True,Theta_Hat)


Number_of_Replications=length(Theta_True);
[Number_of_Examinees, Number_of_Dimensions]=size(Theta_True{1,1});

MSE=zeros(1,Number_of_Dimensions);
Bias=zeros(1,Number_of_Dimensions);
Mean_ED=0;
Correlation_Coefficient=zeros(1,Number_of_Dimensions);

for rep=1:Number_of_Replications
    
    MSE=MSE+(sum((Theta_Hat{rep,1}-Theta_True{rep,1}).^2,1))/Number_of_Examinees;
    Bias=Bias+(sum(Theta_Hat{rep,1}-Theta_True{rep,1},1))/Number_of_Examinees;
    Mean_ED=Mean_ED+mean(sqrt(sum((Theta_Hat{rep,1}-Theta_True{rep,1}).^2,2)));

    Correlation_Coefficient_Temp=zeros(1,Number_of_Dimensions);
    for k=1:(Number_of_Dimensions)
        Correlation_Coefficient_Matrix=corrcoef(Theta_Hat{rep,1}(:,k),Theta_True{rep,1}(:,k));
        Correlation_Coefficient_Temp(1,k)=Correlation_Coefficient_Matrix(1,2);
    end
    Correlation_Coefficient=Correlation_Coefficient+Correlation_Coefficient_Temp;
    
end

MSE=MSE/Number_of_Replications;
Bias=Bias/Number_of_Replications;
Mean_ED=Mean_ED/Number_of_Replications;
Correlation_Coefficient=Correlation_Coefficient/Number_of_Replications;

end
