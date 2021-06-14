function [Sensitivity,Precision,Specificity,Negative,Accuracy] = extractValues(input)

Sensitivity = sum(input(:,1))./( sum(input(:,1))+ sum(input(:,4)))*100;
Precision = sum(input(:,1))./( sum(input(:,1))+ sum(input(:,3)))*100;
Specificity = sum(input(:,2))./( sum(input(:,2))+ sum(input(:,3)))*100;
Negative = sum(input(:,2))./( sum(input(:,2))+ sum(input(:,4)))*100;
Accuracy = (sum(input(:,2))+sum(input(:,1)))./(sum(input(:,2))+ sum(input(:,4))+sum(input(:,1))+ sum(input(:,3)))*100;

if sum(input(:,1))==0,
    Sensitivity = 0;
    Precision   = 0;
elseif sum(input(:,2)) == 0
    Specificity = 0;
    Negative = 0;
end

if sum(input(:,1))==0 && sum(input(:,2))==0
    Accuracy = 0;
end