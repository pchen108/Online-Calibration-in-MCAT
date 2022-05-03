
function [New_Items_Parameter_Estimate_MEM,round]=M_MEM_Online_Calibration(New_Items_Parameter_Estimate_OEM,G_OEM,New_Items_Table,ID_of_Items_Answered,ID_of_New_Items_Answered,V_Matrix,New_Items_V_Matrix,A_Parameter,b_Parameter,Theta_Draw,Prior_Probability,Accuracy_of_Iteration)
% this function is used to calibrate the new items by using M_MEM method

% matrix New_Items_Parameter_Estimate_MEM stores the calibrated results of the new items
% cell matrix G_MEM stores the g(Thetam|Vi) of all new items obtained by M_MEM method
% matrix New_Items_Parameter_Estimate_OEM stores the calibrated results of the new items obtained by M_OEM method
% cell matrix G_OEM stores the g(Thetam|Vi) of all new items obtained by M_OEM method
% cell matrix New_Items_Table records the IDs of the examinees who answered the new items and their responses on the new items
% matrix ID_of_Items_Answered stores the IDs of the operational items which are answered by the examinees
% matrix ID_of_New_Items_Answered stores the IDs of the new items which are answered by the examinees
% matrix V_Matrix stores the response patterns of all examinees on operational items
% matrix New_Items_V_Matrix stores the response patterns of all examinees on new items
% matrix A_Parameter stores all discrimination (a) parameters of all operational items
% column vector b_Parameter stores all b parameters of all operational items
% matrix Theta_Draw stores a random sample drawn from a given distribution
% column vector Prior_Probability stores the prior probability value evaluated at Theta_Draw
% Accuracy_of_Iteration is the precision of iteration we specified when using the Newton-Raphson iterative method and EM algorithm


Number_of_New_Items=length(New_Items_Table);
Number_of_Draws=length(Prior_Probability);

Parameter_MEM=New_Items_Parameter_Estimate_OEM;                % treated as the initial parameter estimates of the new items
Parameter_Hat_MEM=Parameter_MEM;

Number_of_Dimensions=length(Parameter_MEM(1,:));                     % number of parameters need to be estimated, not the number of ability dimensions

G_MEM=G_OEM;                                          

round=1;                                                            % record the number of outer EM cycles (including the first EM cycle)

flag1=1;                                                             % flag1 controls the outer EM cycles
while (flag1==1)
    
    round=round+1;
    disp(['The ',num2str(round),'-th round in M_MEM method!']);
    
    for j=1:Number_of_New_Items                 % visit each new item
    
        Examinee_IDs=(New_Items_Table{j,1}(1,:))';                            % IDs of examinees who answered the current new item
        Number_of_Examinees=length(Examinee_IDs);                         % number of examinees who answered the current new item
    
        % Step 1: E-Step

        % compute likelihood functions (L(Vi|Thetam))
        
        L=zeros(Number_of_Examinees,Number_of_Draws);
    
        for i=1:Number_of_Examinees               % visit each examinee who answered the current new item
            
        	Item_Answered_ID=(ID_of_Items_Answered(Examinee_IDs(i,1),:))';                  % operational items
            V_Answered=(V_Matrix(Examinee_IDs(i,1),:))';
            
            A_Parameter_Answered=A_Parameter(Item_Answered_ID,:);                  % item parameters of the operational items the current examinee answered
            b_Parameter_Answered=b_Parameter(Item_Answered_ID,:);
            
            New_Items_Answered_ID=(ID_of_New_Items_Answered(Examinee_IDs(i,1),:))';          % new items
            New_Items_V_Answered=(New_Items_V_Matrix(Examinee_IDs(i,1),:))';
                        
            New_Items_A_Parameter_Answered=Parameter_MEM(New_Items_Answered_ID,(1:length(A_Parameter(1,:))));        % these two terms are related to Parameter_MEM
            New_Items_b_Parameter_Answered=Parameter_MEM(New_Items_Answered_ID,(length(A_Parameter(1,:))+1));
            
            Total_A_Parameter_Answered=[A_Parameter_Answered;New_Items_A_Parameter_Answered];
            Total_b_Parameter_Answered=[b_Parameter_Answered;New_Items_b_Parameter_Answered];
            Total_V_Answered=[V_Answered;New_Items_V_Answered];
            
            IRFs=1./(1+exp(-Theta_Draw*Total_A_Parameter_Answered').*exp(repmat(Total_b_Parameter_Answered',Number_of_Draws,1)));
        
            L(i,:)=(prod(((IRFs).^(repmat(Total_V_Answered',Number_of_Draws,1))).*((1-IRFs).^(1-(repmat(Total_V_Answered',Number_of_Draws,1)))),2))';
        
        end
        
        % compute the posterior probability of Thetam (mj(Thetam|uij))
    
        % first, compute the posterior probability of Thetam (g(Thetam|Vi))
        
        temp1=L.*G_MEM{j,1};
        g=(temp1)./((1/Number_of_Draws)*repmat(sum(temp1./(repmat(Prior_Probability',Number_of_Examinees,1)),2),1,Number_of_Draws));
        G_MEM{j,1}=g;
                
        % second, compute the likelihood function (Lj(uij|Thetam))
        Response_Pattern=(New_Items_Table{j,1}(2,:))';
        IRFs_New=(Item_Response_Functions_New(Theta_Draw,(Parameter_Hat_MEM(j,:))'))';               % this term is related to Parameter_Hat_MEM
        temp2=repmat(IRFs_New,Number_of_Examinees,1);
        temp3=1-temp2;
        temp4=repmat(Response_Pattern,1,Number_of_Draws);
        Lj=(temp2).^temp4.*(temp3).^(1-temp4);
    
        % third, compute mj
        temp5=(Lj.*G_MEM{j,1})./(repmat(Prior_Probability',Number_of_Examinees,1));
        mj=(temp5)./(repmat(sum(temp5,2),1,Number_of_Draws));
    
        % compute the two artificial data (fjm and rjm)
        fjm=(sum(mj,1))';
        rjm=(sum(temp4.*mj,1))';
        
        % Step 2: M-Step
        
        f=zeros(Number_of_Dimensions,1);
        Df=zeros(Number_of_Dimensions,Number_of_Dimensions);
        
        flag2=1;                                                     % flag2 control the M Steps    
        while (flag2==1)
        
            % first calculate the first order partial derivatives
            Theta_1m=Theta_Draw(:,1);
            Theta_2m=Theta_Draw(:,2);
            Theta_3m=Theta_Draw(:,3);
            IRFs=Item_Response_Functions_New(Theta_Draw,(Parameter_Hat_MEM(j,:))');         % this term is related to Parameter_Hat_MEM
            f(1,1)=sum(Theta_1m.*(rjm-IRFs.*fjm),1);
            f(2,1)=sum(Theta_2m.*(rjm-IRFs.*fjm),1);
            f(3,1)=sum(Theta_3m.*(rjm-IRFs.*fjm),1);
            f(4,1)=-sum(rjm-IRFs.*fjm,1);
                    
            % next calculate the second order partial derivatives
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
                    
            Change_Quantity=Df\f;
            Parameter_Estimate_New=(Parameter_Hat_MEM(j,:))'-Change_Quantity;               % iterative formula
    
            if (max(abs(Change_Quantity))<Accuracy_of_Iteration)
                flag2=0;
            else
                Parameter_Hat_MEM(j,:)=(Parameter_Estimate_New)';
            end

        end
        
        Parameter_Hat_MEM(j,:)=(Parameter_Estimate_New)';
    
    end
    
    Evaluation_Criterion=sum(sum(abs(Parameter_Hat_MEM-Parameter_MEM),1),2)/(Number_of_Dimensions*Number_of_New_Items);
    disp(['The evaluation criterion is ',num2str(Evaluation_Criterion),' !']);
    if (Evaluation_Criterion<Accuracy_of_Iteration)
        flag1=0;
    else
        Parameter_MEM=Parameter_Hat_MEM;
    end

end

New_Items_Parameter_Estimate_MEM=Parameter_Hat_MEM;

end

