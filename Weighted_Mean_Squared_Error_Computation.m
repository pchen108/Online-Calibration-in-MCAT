
function WMSE=Weighted_Mean_Squared_Error_Computation(Item_Parameter_Estimate,Item_Parameter_True,Theta_Draw)



Number_of_Replications=length(Item_Parameter_Estimate);
[Number_of_New_Items, Number_of_Dimensions]=size(Item_Parameter_Estimate{1,1});
Number_of_Draws=length(Theta_Draw{1,1}(:,1));

WMSE=0;
for j=1:Number_of_New_Items
    for r=1:Number_of_Replications
        temp=(Item_Response_Functions_New(Theta_Draw{r,1},(Item_Parameter_True{r,1}(j,:))')-Item_Response_Functions_New(Theta_Draw{r,1},(Item_Parameter_Estimate{r,1}(j,:))')).^2;
        WMSE=WMSE+sum(temp,1)/Number_of_Draws;
    end
end

WMSE=WMSE/(Number_of_New_Items*Number_of_Replications);