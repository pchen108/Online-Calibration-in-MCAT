

%% Matlab simulation program for paper entitled "A Comparative Study of Online Item Calibration in Multidimensional Computerized Adaptive Testing"

clear;
clc;
close all;

% load data
load Condition22;                                             % “N=1800, r=0.5” condition

% define some constants
Sigma_Theta=Sigma_Theta_2;                        % variance-covariance matrix or correlation matrix
Theta_True=True_Theta_22;                           % examinees' ability vectors (100 replications)
Theta_Draw=Theta_Draw_2;                           % 3000 draws used for MC integration
Prior_Probability=Prior_Probability_2;
Random_Matrix=Random_Matrix_2;             % used for Random Design
Score_Matrix=Score_Matrix_22;                    % examinees' responses to the operational items
New_Item_Score_Matrix=Score_Matrix_New_22;              % examinees' responses to the new items

[Number_of_Examinees,Number_of_Dimensions]=size(Theta_True{1,1});                       

Number_of_Operational_Items=length(Score_Matrix{1,1}(1,:));
Number_of_New_Items=length(New_Item_Score_Matrix{1,1}(1,:));
Number_of_Draws=length(Theta_Draw{1,1}(:,1));
Number_of_Replications=length(Theta_True);

MCAT_Test_Length=40;
New_Item_Test_Length=6;
Accuracy_of_Iteration=0.001;
Theta_Initial=[0;0;0];
Number_of_Maximum_Iterations=25;           % used for MLE Estimation Subroutine
Number_of_Quadrature_Nodes=41;              % used for Grid Search method


Theta_Hat=cell(Number_of_Replications,1);                                 % ability vector estimates
Theta_SE=cell(Number_of_Replications,1);
Non_Convergence_for_MLE=cell(Number_of_Replications,1);

New_Items_Method_A_True=cell(Number_of_Replications,1);                              % plugging in the true ability vectors, M-Method A (True)
New_Items_Method_A_Original=cell(Number_of_Replications,1);                        % original M-Method A, M-MethodA (Original)
New_Items_FFMLE_Method_A_Individual=cell(Number_of_Replications,1);     % FFMLE-M-Method A(Individual)
New_Items_FFMLE_Method_A_Mean=cell(Number_of_Replications,1);             % FFMLE-M-Method A(Mean) 
New_Items_M_OEM=cell(Number_of_Replications,1);                                           % M-OEM
New_Items_M_MEM=cell(Number_of_Replications,1);                                          % M-MEM
New_Items_M_OEM_BME=cell(Number_of_Replications,1);                                % M-OEM-BME
New_Items_M_MEM_BME=cell(Number_of_Replications,1);                               % M-MEM-BME

Number_of_EM_Cycles_M_MEM=zeros(Number_of_Replications,1);
Number_of_EM_Cycles_M_MEM_BME=zeros(Number_of_Replications,1);

Elapsed_Time_M_Method_A_True=zeros(Number_of_Replications,1);                 % record the times used for running each method
Elapsed_Time_M_Method_A_Original=zeros(Number_of_Replications,1);
Elapsed_Time_FFMLE_M_Method_A_Individual=zeros(Number_of_Replications,1);
Elapsed_Time_FFMLE_M_Method_A_Mean=zeros(Number_of_Replications,1);
Elapsed_Time_M_OEM=zeros(Number_of_Replications,1);
Elapsed_Time_M_MEM=zeros(Number_of_Replications,1);
Elapsed_Time_M_OEM_BME=zeros(Number_of_Replications,1);
Elapsed_Time_M_MEM_BME=zeros(Number_of_Replications,1);


%% the Entire Process is Repeated /Number_of_Replications/ Times
for rep=1:Number_of_Replications 
    
    Lower_Bound=min(min(Theta_True{rep,1}));                 % lower and upper bound of ability, used for Grid Search method
    Upper_Bound=max(max(Theta_True{rep,1}));
    Theta_All_Possible=Generate_All_Possible_Theta(Number_of_Dimensions,Lower_Bound,Upper_Bound,Number_of_Quadrature_Nodes);

    Mu_AB=Prior_for_New_Items{rep,1}(:,1);
    Sigma_AB=Prior_for_New_Items{rep,1}(:,2:5);
                    
    % Simulate the Entire MCAT Process    
    ID_of_Items_Answered=zeros(Number_of_Examinees,MCAT_Test_Length);         % indicates which operational items are answered by the examinees
    V_Matrix=zeros(Number_of_Examinees,MCAT_Test_Length);                                 % record the response patterns of all examinees
    
    a_Parameter=Item_Parameter_Matrix{rep,1}(:,1:Number_of_Dimensions);              % for convenience
    b_Parameter=Item_Parameter_Matrix{rep,1}(:,Number_of_Dimensions+1);
 
    pointer=0;
    for i=1:Number_of_Examinees                    % visit each examinee
    
        disp(['Rep ',num2str(rep),': the ',num2str(i),'th examinee!']);
        
        Item_Flag=zeros(Number_of_Operational_Items,1);
        U=zeros(MCAT_Test_Length,1);
        V=zeros(MCAT_Test_Length,1);
        Test_Length=0;
        Theta_Estimate=Theta_Initial;
    
        flag=1;
        while (flag==1)
            Item_ID=D_Optimality_Item_Selection_Strategy(Item_Flag,Theta_Estimate,a_Parameter,b_Parameter,Sigma_Theta,U,Test_Length);         % item selection strategy
            Score=Score_Matrix{rep,1}(i,Item_ID);                   % response to item (simulated)
        
            Item_Flag(Item_ID,1)=1;
            Test_Length=Test_Length+1;
            U(Test_Length,1)=Item_ID;
            V(Test_Length,1)=Score;
        
            Theta_Estimate=EAP_Ability_Estimation_Method(a_Parameter,b_Parameter,U,V,Test_Length,Theta_Draw{rep,1});        % ability vector estimation method
        
            if (Test_Length>=MCAT_Test_Length)                   % stopping rule
                flag=0;
            end
        end
    
        [Theta_Estimate,flag1]=MLE_Ability_Estimation_Method(Theta_Estimate,a_Parameter,b_Parameter,U,V,Test_Length,Accuracy_of_Iteration,Number_of_Maximum_Iterations,Lower_Bound,Upper_Bound);
        if (flag1==1)
        	Theta_Estimate=MLE_Grid_Search(Theta_All_Possible,a_Parameter,b_Parameter,U,V,Test_Length);
            pointer=pointer+1;
            Non_Convergence_for_MLE{rep,1}(1,pointer)=i;
        end
            
        Theta_Hat{rep,1}(i,:)=Theta_Estimate';
        Theta_SE{rep,1}(i,:)=(MLE_Standard_Error_Computation(Theta_Estimate,a_Parameter,b_Parameter,U,Test_Length))';

        ID_of_Items_Answered(i,:)=U';
        V_Matrix(i,:)=V';
    
    end
    
    % obtain the Dispersion_Matrix used for FFMLE-M-Method A
    Dispersion_Matrix=cell(Number_of_Examinees,1);

    for i=1:Number_of_Examinees
        Dispersion_Matrix{i,1}=MLE_Sigma_Computation((Theta_Hat{rep,1}(i,:))',a_Parameter,b_Parameter,(ID_of_Items_Answered(i,:))',MCAT_Test_Length);
    end
    
    % Calibrate the New Items by Using Eight Online Calibration Methods (including M-Method A, M-OEM and M-MEM, etc.)
    New_Items_Table=cell(Number_of_New_Items,1);              % record the IDs of the examinees who answered the new item and their responses on the item
    ID_of_New_Items_Answered=zeros(Number_of_Examinees,New_Item_Test_Length);            % indicate which new items the examinees answered
    New_Item_V_Matrix=zeros(Number_of_Examinees,New_Item_Test_Length);
    
    for j=1:Number_of_New_Items
        k=0;
        for i=1:Number_of_Examinees
            if (Random_Matrix{rep,1}(i,j)==1)              
                k=k+1;
                New_Items_Table{j,1}(1,k)=i;                            % record which examinees answered this new item
                New_Items_Table{j,1}(2,k)=New_Item_Score_Matrix{rep,1}(i,j);                 % item responses to this new item
            end
        end
    end
    
    for i=1:Number_of_Examinees
        k=0;
        for j=1:Number_of_New_Items
            if (Random_Matrix{rep,1}(i,j)==1)            
                k=k+1;
                ID_of_New_Items_Answered(i,k)=j;
            end
        end
    end

    for i=1:Number_of_Examinees
    	New_Item_V_Matrix(i,:)=New_Item_Score_Matrix{rep,1}(i,ID_of_New_Items_Answered(i,:));
    end

    Parameter_Initial=zeros(Number_of_New_Items,(Number_of_Dimensions+1));
    Success_Rate=(sum(New_Item_Score_Matrix{rep,1},1)/Number_of_Examinees)';
    Parameter_Initial(:,(Number_of_Dimensions+1))=-norminv(Success_Rate,0,1);
    Parameter_Initial(:,1:Number_of_Dimensions)=0.5*ones(Number_of_New_Items,Number_of_Dimensions);           % initial values for a is 0.5
        
    disp('Start to Run M-Method A(True) and M-Method A(Original)!');
    tic;
    New_Items_Method_A_True{rep,1}=M_Method_A_Online_Calibration(Theta_True{rep,1},New_Items_Table,Parameter_Initial,Accuracy_of_Iteration);
    Elapsed_Time_M_Method_A_True(rep,1)=toc;
    tic;
    New_Items_Method_A_Original{rep,1}=M_Method_A_Online_Calibration(Theta_Hat{rep,1},New_Items_Table,Parameter_Initial,Accuracy_of_Iteration);
    Elapsed_Time_M_Method_A_Original(rep,1)=toc;
    
    disp('Start to Run FFMLE-M-Method A(Individual) and FFMLE-M-Method A(Mean)!');
    tic;
    New_Items_FFMLE_Method_A_Individual{rep,1}=FFMLE_M_Method_A_Online_Calibration_Individual(Theta_Hat{rep,1},New_Items_Table,New_Items_Method_A_Original{rep,1},Dispersion_Matrix,Parameter_Initial,Accuracy_of_Iteration);
    Elapsed_Time_FFMLE_M_Method_A_Individual(rep,1)=toc;
    tic;
    New_Items_FFMLE_Method_A_Mean{rep,1}=FFMLE_M_Method_A_Online_Calibration_Mean(Theta_Hat{rep,1},New_Items_Table,New_Items_Method_A_Original{rep,1},Dispersion_Matrix,Parameter_Initial,Accuracy_of_Iteration);
    Elapsed_Time_FFMLE_M_Method_A_Mean(rep,1)=toc;
    
    disp('Start to Run M-OEM and M-MEM!');
    tic;
    [New_Items_M_OEM{rep,1},G_OEM]=M_OEM_Online_Calibration(New_Items_Table,ID_of_Items_Answered,V_Matrix,a_Parameter,b_Parameter,Parameter_Initial,Theta_Draw{rep,1},Prior_Probability{rep,1},Accuracy_of_Iteration);
    Elapsed_Time_M_OEM(rep,1)=toc;
    tic;
    [New_Items_M_MEM{rep,1},Number_of_EM_Cycles_M_MEM(rep,1)]=M_MEM_Online_Calibration(New_Items_M_OEM{rep,1},G_OEM,New_Items_Table,ID_of_Items_Answered,ID_of_New_Items_Answered,V_Matrix,New_Item_V_Matrix,a_Parameter,b_Parameter,Theta_Draw{rep,1},Prior_Probability{rep,1},Accuracy_of_Iteration);
    Elapsed_Time_M_MEM(rep,1)=toc;
    
    disp('Start to Run M-OEM-BME and M-MEM-BME!');
    tic;
    [New_Items_M_OEM_BME{rep,1},G_OEM_BME]=M_OEM_BME_Online_Calibration(New_Items_Table,ID_of_Items_Answered,V_Matrix,a_Parameter,b_Parameter,Parameter_Initial,Theta_Draw{rep,1},Prior_Probability{rep,1},Accuracy_of_Iteration,Mu_AB,Sigma_AB);
    Elapsed_Time_M_OEM_BME(rep,1)=toc;
    tic;
    [New_Items_M_MEM_BME{rep,1},Number_of_EM_Cycles_M_MEM_BME(rep,1)]=M_MEM_BME_Online_Calibration(New_Items_M_OEM_BME{rep,1},G_OEM_BME,New_Items_Table,ID_of_Items_Answered,ID_of_New_Items_Answered,V_Matrix,New_Item_V_Matrix,a_Parameter,b_Parameter,Theta_Draw{rep,1},Prior_Probability{rep,1},Accuracy_of_Iteration,Mu_AB,Sigma_AB);
    Elapsed_Time_M_MEM_BME(rep,1)=toc;
    
end


%% Compute the Evaluation Criteria regarding Ability Vector Estimation
[MSE_Theta,Bias_Theta,Correlation_Coefficient_Theta,Mean_ED]=Evaluation_Criterion_About_Ability_Parameter(Theta_True,Theta_Hat);
Results_Theta=[MSE_Theta,Bias_Theta,Correlation_Coefficient_Theta,Mean_ED];
    
% Compute the Evaluation Criteria regarding Item Parameter Estimation
% M-Method A (True)
[MSE_True,Bias_True,Correlation_Coefficient_True]=Evaluation_Criterion_About_Item_Parameter(New_Items_Method_A_True,New_Item_Parameter_Matrix);
WMSE_True=Weighted_Mean_Squared_Error_Computation(New_Items_Method_A_True,New_Item_Parameter_Matrix,Theta_Draw);
Results_Hecheng_True=[MSE_True,Bias_True,Correlation_Coefficient_True,WMSE_True];

% M-Method A (Original)
[MSE_Original,Bias_Original,Correlation_Coefficient_Original]=Evaluation_Criterion_About_Item_Parameter(New_Items_Method_A_Original,New_Item_Parameter_Matrix);
WMSE_Original=Weighted_Mean_Squared_Error_Computation(New_Items_Method_A_Original,New_Item_Parameter_Matrix,Theta_Draw);
Results_Original=[MSE_Original,Bias_Original,Correlation_Coefficient_Original,WMSE_Original];

% FFMLE-M-Method A (Individual)
[MSE_FFMLE_Individual,Bias_FFMLE_Individual,Correlation_Coefficient_FFMLE_Individual]=Evaluation_Criterion_About_Item_Parameter(New_Items_FFMLE_Method_A_Individual,New_Item_Parameter_Matrix);
WMSE_FFMLE_Individual=Weighted_Mean_Squared_Error_Computation(New_Items_FFMLE_Method_A_Individual,New_Item_Parameter_Matrix,Theta_Draw);
Results_FFMLE_Individual=[MSE_FFMLE_Individual,Bias_FFMLE_Individual,Correlation_Coefficient_FFMLE_Individual,WMSE_FFMLE_Individual];

% FFMLE-M-Method A (Mean)
[MSE_FFMLE_Mean,Bias_FFMLE_Mean,Correlation_Coefficient_FFMLE_Mean]=Evaluation_Criterion_About_Item_Parameter(New_Items_FFMLE_Method_A_Mean,New_Item_Parameter_Matrix);
WMSE_FFMLE_Mean=Weighted_Mean_Squared_Error_Computation(New_Items_FFMLE_Method_A_Mean,New_Item_Parameter_Matrix,Theta_Draw);
Results_FFMLE_Mean=[MSE_FFMLE_Mean,Bias_FFMLE_Mean,Correlation_Coefficient_FFMLE_Mean,WMSE_FFMLE_Mean];

% M-OEM
[MSE_OEM,Bias_OEM,Correlation_Coefficient_OEM]=Evaluation_Criterion_About_Item_Parameter(New_Items_M_OEM,New_Item_Parameter_Matrix);
WMSE_OEM=Weighted_Mean_Squared_Error_Computation(New_Items_M_OEM,New_Item_Parameter_Matrix,Theta_Draw);
Results_OEM=[MSE_OEM,Bias_OEM,Correlation_Coefficient_OEM,WMSE_OEM];
    
% M-MEM
[MSE_MEM,Bias_MEM,Correlation_Coefficient_MEM]=Evaluation_Criterion_About_Item_Parameter(New_Items_M_MEM,New_Item_Parameter_Matrix);
WMSE_MEM=Weighted_Mean_Squared_Error_Computation(New_Items_M_MEM,New_Item_Parameter_Matrix,Theta_Draw);
Results_MEM=[MSE_MEM,Bias_MEM,Correlation_Coefficient_MEM,WMSE_MEM];
    
% M-OEM-BME
[MSE_OEM_BME,Bias_OEM_BME,Correlation_Coefficient_OEM_BME]=Evaluation_Criterion_About_Item_Parameter(New_Items_M_OEM_BME,New_Item_Parameter_Matrix);
WMSE_OEM_BME=Weighted_Mean_Squared_Error_Computation(New_Items_M_OEM_BME,New_Item_Parameter_Matrix,Theta_Draw);
Results_OEM_BME=[MSE_OEM_BME,Bias_OEM_BME,Correlation_Coefficient_OEM_BME,WMSE_OEM_BME];
    
% M-MEM-BME
[MSE_MEM_BME,Bias_MEM_BME,Correlation_Coefficient_MEM_BME]=Evaluation_Criterion_About_Item_Parameter(New_Items_M_MEM_BME,New_Item_Parameter_Matrix);
WMSE_MEM_BME=Weighted_Mean_Squared_Error_Computation(New_Items_M_MEM_BME,New_Item_Parameter_Matrix,Theta_Draw);
Results_MEM_BME=[MSE_MEM_BME,Bias_MEM_BME,Correlation_Coefficient_MEM_BME,WMSE_MEM_BME];

% record the average times used for running each online calibration method
Time_M_Method_A_True=mean(Elapsed_Time_M_Method_A_True);
Time_M_Method_A_Original=mean(Elapsed_Time_M_Method_A_Original);
Time_FFMLE_M_Method_A_Individual=mean(Elapsed_Time_FFMLE_M_Method_A_Individual);
Time_FFMLE_M_Method_A_Mean=mean(Elapsed_Time_FFMLE_M_Method_A_Mean);
Time_M_OEM=mean(Elapsed_Time_M_OEM);
Time_M_MEM=mean(Elapsed_Time_M_MEM);
Time_M_OEM_BME=mean(Elapsed_Time_M_OEM_BME);
Time_M_MEM_BME=mean(Elapsed_Time_M_MEM_BME);

% record the average number of EM cycles used for running M-MEM and M-MEM-BME method
Average_Number_of_EM_Cycles_M_MEM=mean(Number_of_EM_Cycles_M_MEM);
Average_Number_of_EM_Cycles_M_MEM_BME=mean(Number_of_EM_Cycles_M_MEM_BME);


%% Store Some Results
save Evaluation_Results.mat Results_Theta Results_True Results_Original Results_FFMLE_Individual Results_FFMLE_Mean Results_OEM Results_MEM Results_OEM_BME Results_MEM_BME
save Number_of_EM_Cycles.mat Number_of_EM_Cycles_M_MEM Number_of_EM_Cycles_M_MEM_BME Average_Number_of_EM_Cycles_M_MEM Average_Number_of_EM_Cycles_M_MEM_BME
save Time_Used.mat Time_M_Method_A_True Time_M_Method_A_Original Time_FFMLE_M_Method_A_Individual Time_FFMLE_M_Method_A_Mean Time_M_OEM Time_M_MEM Time_M_OEM_BME Time_M_MEM_BME
save Ability_Related.mat Theta_Hat Theta_SE
save Calibration_Results.mat New_Items_Method_A_True New_Items_Method_A_Original New_Items_FFMLE_Method_A_Individual New_Items_FFMLE_Method_A_Mean New_Items_M_OEM New_Items_M_MEM New_Items_M_OEM_BME New_Items_M_MEM_BME
save Non_Convergence_for_MLE.mat Non_Convergence_for_MLE


%%
load Non_Convergence_for_MLE.mat;
load Ability_Related.mat;
load Condition22.mat;                                     % (CHANGE)
Theta_True=True_Theta_22;                           % (CHANGE)
[MSE_Theta_1,Bias_Theta_1,Correlation_Coefficient_Theta_1,Mean_ED_1]=Evaluation_Criterion_About_Ability_Parameter_1(Theta_True,Theta_Hat,Non_Convergence_for_MLE);
Results_Hecheng_Theta_1=[MSE_Theta_1,Bias_Theta_1,Correlation_Coefficient_Theta_1,Mean_ED_1];
