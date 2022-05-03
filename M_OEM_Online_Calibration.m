
function [New_Items_Parameter_Estimate,G_OEM]=M_OEM_Online_Calibration(New_Items_Table,ID_of_Items_Answered,V_Matrix,A_Parameter,b_Parameter,Parameter_Initial,Theta_Draw,Prior_Probability,Accuracy_of_Iteration)
% this function is used to calibrate the new items by using M_OEM method

% matrix New_Items_Parameter_Estimate stores the estimated parameters of the new items
% cell matrix G_OEM stores the g(Thetam|Vi) obtained by M_OEM method 
% cell Matrix New_Items_Table records the IDs of the examinees who answered the new items and their responses on the items
% matrix ID_of_Items_Answered stores the IDs of the items which are answered by the examinees
% matrix V_Matrix stores the response patterns of all examinees
% matrix A_Parameter stores all discrimination parameters of all items
% column vector b_Parameter stores all b parameters of all items
% matrix Parameter_Initial stores the initial parameter values 
% matrix Theta_Draw stores a random sample drawn from a given distribution
% column vector Prior_Probability stores the prior probability value evaluated at Theta_Draw
% Accuracy_of_Iteration is the precision of iteration we specified when using the Newton-Raphson iterative method


Number_of_New_Items=length(New_Items_Table);
Number_of_Draws=length(Theta_Draw(:,1));

Number_of_Dimensions=length(A_Parameter(1,:))+1;                        % number of parameters need to be estimated, rather than the number of ability dimensions

New_Items_Parameter_Estimate=zeros(Number_of_New_Items,Number_of_Dimensions);

% G_OEM=zeros(Number_of_Draws,Number_of_New_Items);       % row represents each draw, column represents each new item
G_OEM=cell(Number_of_New_Items,1);

for j=1:Number_of_New_Items
    
    Parameter_Estimate=(Parameter_Initial(j,:))';             % initial parameter values of each new item

    Examinee_IDs=(New_Items_Table{j,1}(1,:))';           % IDs of examinees who answered the current new item
    Number_of_Examinees=length(Examinee_IDs);        % number of examinees who answered the current new item
    
    % Step 1: E-Step

    % compute likelihood functions (L(Vi|Thetam))
    L=zeros(Number_of_Examinees,Number_of_Draws);
    
    for i=1:Number_of_Examinees                   % visit each examinee who answered the current new item
        
        Item_Answered_ID=(ID_of_Items_Answered(Examinee_IDs(i,1),:))';
        V_Answered=(V_Matrix(Examinee_IDs(i,1),:))';
        
        A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);             % item parameters of the operational items the current examinee answered
        b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);
        
        IRFs=1./(1+exp(-Theta_Draw*A_Parameter_Answered').*exp(repmat(b_Parameter_Answered',Number_of_Draws,1)));
        
        L(i,:)=(prod(((IRFs).^(repmat(V_Answered',Number_of_Draws,1))).*((1-IRFs).^(1-(repmat(V_Answered',Number_of_Draws,1)))),2))';
        
    end
        
    % compute the posterior probability of Thetam (mj(Thetam|uij))
    
    % first, compute the posterior probability of Thetam (g(Thetam|Vi))
    temp1=L.*(repmat(Prior_Probability',Number_of_Examinees,1));
    g=(temp1)./((1/Number_of_Draws)*repmat(sum(L,2),1,Number_of_Draws));
    G_OEM{j,1}=g;
    
    % second, compute the likelihood function (Lj(uij|Thetam))
    Response_Pattern=(New_Items_Table{j,1}(2,:))';
    IRFs_New=(1./(1+exp(-Theta_Draw*Parameter_Estimate(1:(Number_of_Dimensions-1),1)).*exp(repmat(Parameter_Estimate(Number_of_Dimensions,1),Number_of_Draws,1))))';
    temp2=repmat(IRFs_New,Number_of_Examinees,1);
    temp3=1-temp2;
    temp4=repmat(Response_Pattern,1,Number_of_Draws);
    Lj=(temp2).^temp4.*(temp3).^(1-temp4);
    
    % third, compute mj
    temp5=(Lj.*g)./(repmat(Prior_Probability',Number_of_Examinees,1));
    mj=(temp5)./(repmat(sum(temp5,2),1,Number_of_Draws));
        
    % next, compute the two artificial data (fjm and rjm)
    fjm=(sum(mj,1))';
    rjm=(sum(temp4.*mj,1))';
        
    % Step 2: M-Step
        
    f=zeros(Number_of_Dimensions,1);
    Df=zeros(Number_of_Dimensions,Number_of_Dimensions);
    flag=1;
    
    while (flag==1)
        
        % first, calculate the first order partial derivatives
        Theta_1m=Theta_Draw(:,1);
        Theta_2m=Theta_Draw(:,2);
        Theta_3m=Theta_Draw(:,3);
        IRFs=1./(1+exp(-Theta_Draw*Parameter_Estimate(1:(Number_of_Dimensions-1),1)).*exp(repmat(Parameter_Estimate(Number_of_Dimensions,1),Number_of_Draws,1)));
        f(1,1)=sum(Theta_1m.*(rjm-IRFs.*fjm),1);
        f(2,1)=sum(Theta_2m.*(rjm-IRFs.*fjm),1);
        f(3,1)=sum(Theta_3m.*(rjm-IRFs.*fjm),1);
        f(4,1)=-sum(rjm-IRFs.*fjm,1);
                
        % second, calculate the second order partial derivatives
        Df(1,1)=-sum((Theta_1m).^2.*fjm.*IRFs.*(1-IRFs),1);
        Df(2,2)=-sum((Theta_2m).^2.*fjm.*IRFs.*(1-IRFs),1);
        Df(3,3)=-sum((Theta_3m).^2.*fjm.*IRFs.*(1-IRFs),1);
        Df(4,4)=-sum(fjm.*IRFs.*(1-IRFs),1);
        Df(1,2)=-sum(Theta_1m.*Theta_2m.*fjm.*IRFs.*(1-IRFs),1);
        Df(1,3)=-sum(Theta_1m.*Theta_3m.*fjm.*IRFs.*(1-IRFs),1);
        Df(1,4)=sum(Theta_1m.*fjm.*IRFs.*(1-IRFs),1);
        Df(2,3)=-sum(Theta_2m.*Theta_3m.*fjm.*IRFs.*(1-IRFs),1);
        Df(2,4)=sum(Theta_2m.*fjm.*IRFs.*(1-IRFs),1);
        Df(3,4)=sum(Theta_3m.*fjm.*IRFs.*(1-IRFs),1);
        Df(2,1)=Df(1,2);
        Df(3,1)=Df(1,3);
        Df(3,2)=Df(2,3);
        Df(4,1)=Df(1,4);
        Df(4,2)=Df(2,4);
        Df(4,3)=Df(3,4);
                
        % next, Newton-Raphson iteration
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

