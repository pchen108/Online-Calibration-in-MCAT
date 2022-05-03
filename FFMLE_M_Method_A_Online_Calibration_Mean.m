
function New_Items_Parameter_Estimate=FFMLE_M_Method_A_Online_Calibration_Mean(Theta_Estimate,New_Items_Table,Beta_Hat,Dispersion_Matrix,Parameter_Initial,Accuracy_of_Iteration)
% this function is used to calibrate the new items by using M_Method_A method

% matrix New_Items_Parameter_Estimate stores the estimated parameters of the new items
% matrix Theta_Estimate stores the estimated ability vectors of all examinees which are considered as the true ability vector in M_Method_A method
% cell Matrix New_Items_Table records the IDs of the examinees who answered the new items and their responses on the items
% matrix Parameter_Initial stores the initial parameter values 
% Accuracy_of_Iteration is the precision of iteration we specified when using the Newton-Raphson iterative method


[Number_of_New_Items,Number_of_Dimensions]=size(Parameter_Initial);                  % number of parameters need to be estimated, rather than the number of ability dimensions

New_Items_Parameter_Estimate=zeros(Number_of_New_Items,Number_of_Dimensions);

for j=1:Number_of_New_Items
    
    Parameter_Estimate=(Parameter_Initial(j,:))';     % initial parameter values

    Examinee_IDs=(New_Items_Table{j,1}(1,:))';
    Response_Pattern=(New_Items_Table{j,1}(2,:))';
    Theta_True_f=zeros(length(Examinee_IDs),Number_of_Dimensions-1);
    Dispersion_Matrix_Mean=zeros(Number_of_Dimensions-1,Number_of_Dimensions-1);
    Beta_j=(Beta_Hat(j,:))';
    
    for i=1:(length(Examinee_IDs))
        Dispersion_Matrix_Mean=Dispersion_Matrix_Mean+Dispersion_Matrix{Examinee_IDs(i,1),1};
    end
    Dispersion_Matrix_Mean=Dispersion_Matrix_Mean/(length(Examinee_IDs));
    
    for i=1:(length(Examinee_IDs))
        ID=Examinee_IDs(i,1);
        Theta_True_f(i,:)=((Theta_Estimate(ID,:))'+(Response_Pattern(i,1)-Item_Response_Functions_New(Theta_Estimate(ID,:),Beta_j))*Dispersion_Matrix_Mean*Beta_j(1:(Number_of_Dimensions-1),1))';
    end
    
    f=zeros(Number_of_Dimensions,1);
    Df=zeros(Number_of_Dimensions,Number_of_Dimensions);
    flag=1;
    
    while (flag==1)
        
        % first calculate the first order partial derivatives
        IRFs_New=1./(1+exp(-Theta_True_f*Parameter_Estimate(1:(Number_of_Dimensions-1),1)).*exp(repmat(Parameter_Estimate(Number_of_Dimensions,1),length(Theta_True_f(:,1)),1)));
        f(1,1)=sum(Theta_True_f(:,1).*(Response_Pattern-IRFs_New),1);
        f(2,1)=sum(Theta_True_f(:,2).*(Response_Pattern-IRFs_New),1);
        f(3,1)=sum(Theta_True_f(:,3).*(Response_Pattern-IRFs_New),1);
        f(4,1)=-sum((Response_Pattern-IRFs_New),1);
        
        % next calculate the second order partial derivatives
        Df(1,1)=-sum(Theta_True_f(:,1).^2.*IRFs_New.*(1-IRFs_New),1);
        Df(2,2)=-sum(Theta_True_f(:,2).^2.*IRFs_New.*(1-IRFs_New),1);
        Df(3,3)=-sum(Theta_True_f(:,3).^2.*IRFs_New.*(1-IRFs_New),1);
        Df(4,4)=-sum(IRFs_New.*(1-IRFs_New),1);
        Df(1,2)=-sum(Theta_True_f(:,1).*Theta_True_f(:,2).*IRFs_New.*(1-IRFs_New),1);
        Df(1,3)=-sum(Theta_True_f(:,1).*Theta_True_f(:,3).*IRFs_New.*(1-IRFs_New),1);
        Df(1,4)=sum(Theta_True_f(:,1).*IRFs_New.*(1-IRFs_New),1);
        Df(2,3)=-sum(Theta_True_f(:,2).*Theta_True_f(:,3).*IRFs_New.*(1-IRFs_New),1);
        Df(2,4)=sum(Theta_True_f(:,2).*IRFs_New.*(1-IRFs_New),1);
        Df(3,4)=sum(Theta_True_f(:,3).*IRFs_New.*(1-IRFs_New),1);
        Df(2,1)=Df(1,2);
        Df(3,1)=Df(1,3);
        Df(3,2)=Df(2,3);
        Df(4,1)=Df(1,4);
        Df(4,2)=Df(2,4);
        Df(4,3)=Df(3,4);
        
        Change_Quantity=Df\f;
        Parameter_Estimate_New=Parameter_Estimate-Change_Quantity;               % iterative formula
    
        if (max(abs(Change_Quantity))<Accuracy_of_Iteration)
            flag=0;
        else
            Parameter_Estimate=Parameter_Estimate_New;
        end

    end
    
    New_Items_Parameter_Estimate(j,:)=Parameter_Estimate_New';
    
end

end

