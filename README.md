# Online-Calibration-in-MCAT

Matlab simulation program for the 2017 JEBS paper entitled "A Comparative Study of Online Item Calibration in Multidimensional Computerized Adaptive Testing", including one main function (Main_Procedure.m), 19 subroutines (e.g., M_Method_A_Online_Calibration.m), and one dataset file (Condition 22.mat).

Note that: the dataset file is larger than 25MB and can not be uploaded. Please contact me if you need the dataset.

The functions of the 19 subroutines are briefly described as follows:

1. D_Optimality_Item_Selection_Strategy.m: select the next item using the D-Optimality item selection strategy
2. EAP_Ability_Estimation_Method.m: estimate the examinees' ability vectors using EAP method
3. Evaluation_Criterion_About_Ability_Parameter.m: compute the evaluation index value for theta
4. Evaluation_Criterion_About_Ability_Parameter_1.m: compute the evaluation index value for theta considering Non_Convergence case for MLE
5. Evaluation_Criterion_About_Item_Parameter.m: compute the evaluation index value for item parameter
6. FFMLE_M_Method_A_Online_Calibration_Individual.m: calibrate the new items using FFMLE-M-Method A(Individual)
7. FFMLE_M_Method_A_Online_Calibration_Mean.m: calibrate the new items using FFMLE-M-Method A(Mean)
8. Generate_All_Possible_Theta.m: generate all possible ability vectors for the MLE_Grid_Search subroutine
9. Item_Response_Functions_New.m: compute the IRFs of the current new item according to M2PLM
10. M_MEM_BME_Online_Calibration.m: calibrate the new items via M-MEM-BME
11. M_MEM_Online_Calibration.m: calibrate the new items via M-MEM
12. M_Method_A_Online_Calibration.m: calibrate the new items via M-Method A
13. M_OEM_BME_Online_Calibration.m: calibrate the new items via M-OEM-BME
14. M_OEM_Online_Calibration.m: calibrate the new items via M-OEM
15. MLE_Ability_Estimation_Method.m: estimate examinees' ability vectors using MLE method
16. MLE_Grid_Search.m: find the MLE estimate using grid search method
17. MLE_Sigma_Computation.m: calculate the dispersion matrix used for FFMLE-M-Method A method
18. MLE_Standard_Error_Computation.m: calculate the SE for MLE method
19. Weighted_Mean_Squared_Error_Computation.m: calculate the WMSE index
