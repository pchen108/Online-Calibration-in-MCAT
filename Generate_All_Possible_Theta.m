
function Theta_All_Possible=Generate_All_Possible_Theta(Number_of_Dimensions,Lower_Bound,Upper_Bound,Number_of_Quadrature_Nodes)
% this function is used to generate all possible ability vectors for the MLE_Grid_Search ability estimation

% matrix Theta_All_Possible stores all possible ability vectors
% Number_of_Dimensions is the number of dimensions
% Lower_Bound is the lower bound of ability value we set when taking quadrature nodes
% Upper_Bound is the upper bound of ability value we set when taking quadrature nodes


Number_of_Quadrature_Nodes_Total=Number_of_Quadrature_Nodes^(Number_of_Dimensions);

Step_Size=(Upper_Bound-Lower_Bound)/(Number_of_Quadrature_Nodes-1);        % obtain the step size

Theta_All_Possible=zeros(Number_of_Dimensions,Number_of_Quadrature_Nodes_Total);           % preallocating for speed

% compute all possible ability vectors
for i=1:Number_of_Quadrature_Nodes
    value1=Lower_Bound+Step_Size*(i-1);
    for j=1:Number_of_Quadrature_Nodes
        value2=Lower_Bound+Step_Size*(j-1);
        for k=1:Number_of_Quadrature_Nodes
            value3=Lower_Bound+Step_Size*(k-1);
            value=[value1;value2;value3];
            Theta_All_Possible(:,(Number_of_Quadrature_Nodes)^2*(i-1)+(Number_of_Quadrature_Nodes)*(j-1)+k)=value;
        end
    end
end

end
