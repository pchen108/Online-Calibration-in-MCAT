
function [MLE_Estimate,flag1]=MLE_Ability_Estimation_Method(Theta_Estimate,A_Parameter,b_Parameter,U,V,Test_Length,Accuracy_of_Iteration,Number_of_Maximum_Iterations,Lower_Bound,Upper_Bound)
% this function is used to estimate examinees' ability vectors using Maximum Likelihood Estimation (MLE) method

% column vector MLE_Estimate stores the MLE estimate of examinee
% column vector Theta_Estimate stores the estimated ability vector of examinee
% matrix A_Parameter stores the discrimination parameters of all items
% column vector b_Parameter stores the b parameters of all items
% column vector U stores the IDs of items which the current examinee has answered
% column vector V records the responses to the items which the current examinee has answered
% Test_Length is current test length
% Accuracy_of_Iteration is the precision of iteration we specified when using the Newton-Raphson iterative method



Number_of_Dimensions=length(A_Parameter(1,:));

Item_Answered_ID=U(1:Test_Length,:);
A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);
b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);

flag=1;
flag1=0;
Number_of_Iterations=0;

while (flag==1)

    Number_of_Iterations=Number_of_Iterations+1;
    if (Number_of_Iterations>Number_of_Maximum_Iterations)           % the maximum number of iterations
        flag1=1;                                                     % abnormally Exit 1
        break;
    end
    
    % compute the item response functions of the Test_Length items at Theta_Estimate
    IRFs=1./(1+exp(-A_Parameter_Answered*Theta_Estimate).*exp(b_Parameter_Answered));
    
    % Calculate the First Order Partial Derivatives (column vector of Number_of_Dimensions*1)
    Response_Pattern=V(1:Test_Length,:);
    temp1=Response_Pattern-IRFs;
    temp2=repmat(temp1,1,Number_of_Dimensions);
    temp3=(sum(A_Parameter_Answered.*temp2,1))';

    First_Order_Partial_Derivatives=temp3;

    % Calculate the Second Order Partial Derivatives
    Second_Order_Partial_Derivatives=zeros(Number_of_Dimensions,Number_of_Dimensions);

    % first compute the main-diagonal elements
    temp4=IRFs.*(1-IRFs);
    temp5=repmat(temp4,1,Number_of_Dimensions);
    temp6=(A_Parameter_Answered).^2;
    temp7=(-sum(temp6.*temp5,1))';

    for k=1:Number_of_Dimensions
        Second_Order_Partial_Derivatives(k,k)=temp7(k,1);
    end

    % then compute the off-diagonal elements
    for k=1:Number_of_Dimensions
        for l=1:Number_of_Dimensions
            if (k<l)
                Second_Order_Partial_Derivatives(k,l)=-sum(A_Parameter_Answered(:,k).*A_Parameter_Answered(:,l).*temp4,1);
                Second_Order_Partial_Derivatives(l,k)=Second_Order_Partial_Derivatives(k,l);
            end
        end
    end

    Condition_Number=cond(Second_Order_Partial_Derivatives);
    if (Condition_Number>75)
        u=0;
        for k=1:Number_of_Dimensions
            for l=(k+1):Number_of_Dimensions
                u=u+abs(Second_Order_Partial_Derivatives(k,l));
            end
        end
        Second_Order_Partial_Derivatives=Second_Order_Partial_Derivatives+u*eye(Number_of_Dimensions);
    end
    
    Change_Quantity=Second_Order_Partial_Derivatives\First_Order_Partial_Derivatives;
    Theta_Estimate_New=Theta_Estimate-Change_Quantity;               % iterative formula
    
    for k=1:Number_of_Dimensions
        if (Theta_Estimate_New(k,1)<Lower_Bound)
            Theta_Estimate_New(k,1)=Lower_Bound;
        elseif (Theta_Estimate_New(k,1)>Upper_Bound)
            Theta_Estimate_New(k,1)=Upper_Bound;
        end
    end
    
    if (max(abs(Change_Quantity))<Accuracy_of_Iteration)
        flag=0;
    else
        Theta_Estimate=Theta_Estimate_New;
    end
    
end

MLE_Estimate=Theta_Estimate_New;

end
