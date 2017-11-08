%HW6

% Problem 1. Curve fitting. 
% Part 1. Take 10 x values on the interval 0 < x <= 10 and then create y
% values from the x values by plugging the x values into a third order
% polynomial of your choice. Add random noise to the data by choosing a random number
% in the interval [-D, D]. Start with D = 1. Plot your data in the x-y plane.
xx = 10*rand(10,1);
xx = sort(xx);


yy = (xx.^3 + xx +1) + 2*rand(10,1)-1;
% my equation is (x.^3 + x +1), 2*rand(10,1)-1 adding noise 
scatter(xx,yy); 

% Part 2. Fit your data with polynomials from order 1 to 9. Plot the fitted
% polynomials together with the data. 
for order =1:9
    order
    [coeff, s] = polyfit(xx, yy, order);
    figure (order)
    
    plot(xx,yy, 'r.', 'MarkerSize', 24); hold on; 
    plot (xx, polyval(coeff,xx), 'k-', 'LineWidth',3);
    title (order)
    hold off; 
end 
    

% Part 3. On a separate plot, plot the R^2 and adjusted R^2 as a function
% of the order of the polynomial. 
figure(1)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2= fit_metric.rsquare
    
    hold on 
    plot(order, r2, 'r.')
    
    title("r2");
    
    
    
end 
hold off 

figure(2)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2adj = fit_metric.adjrsquare
    hold on 
    plot(order, r2adj, 'b.')
    title("r2adj")
    
    
end
hold off 

% Part 4. Repeat parts 1 - 3 for D = 10 and D = 1000. Comment on the
% results. 

%%D=10
xx = 10*rand(10,1);
xx = sort(xx);
yy = (xx.^3 + xx +1) + 20*rand(10,1)-10;
% my equation is (x.^3 + x +1), 2*rand(10,1)-1 adding noise 
scatter(xx,yy); 
%Fitting Data
for order =1:9
    order
    [coeff, s] = polyfit(xx, yy, order);
    figure (order)
    plot(xx,yy, 'r.', 'MarkerSize', 24); hold on; 
    plot (xx, polyval(coeff,xx), 'k-', 'LineWidth',3);
    title (order)
    hold off; 
end 

%checking fitting
figure(19)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2= fit_metric.rsquare
    
    hold on 
    plot(order, r2, 'r.')
    
    title("r2");
    
end 
hold off 

figure(20)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2adj = fit_metric.adjrsquare
    hold on 
    plot(order, r2adj, 'b.')
    title("r2adj")
    
end
hold off 

%%D=1000
xx = 10*rand(10,1);
xx = sort(xx);
yy = (xx.^3 + xx +1) + 2000*rand(10,1)-1000;
% my equation is (x.^3 + x +1), 2*rand(10,1)-1 adding noise 
scatter(xx,yy); 
%Fitting Data
for order =1:9
    order
    [coeff, s] = polyfit(xx, yy, order);
    figure (order)
    plot(xx,yy, 'r.', 'MarkerSize', 24); hold on; 
    plot (xx, polyval(coeff,xx), 'k-', 'LineWidth',3);
    title (order)
    hold off; 
end 

%checking fitting
figure(19)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2= fit_metric.rsquare
    
    hold on 
    plot(order, r2, 'r.')
    
    title("r2");
    
end 
hold off 

figure(20)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2adj = fit_metric.adjrsquare
    hold on 
    plot(order, r2adj, 'b.')
    title("r2adj")
    
end
hold off 
% Comments: When the noise is big, the model cannot explain the variation
% as seen in the decreased r^2 values. Having a higher order model can fit
% the higher variation better. 

% Part 5. Now repeat parts 1-3 but take 100 x values on the interval 0 < x <=
% 10. Comment on the results. 

%%D=10
xx = 10*rand(100,1);
xx = sort(xx);
yy = (xx.^3 + xx +1) + 20*rand(100,1)-10;
% my equation is (x.^3 + x +1), 2*rand(10,1)-1 adding noise 
scatter(xx,yy); 
%Fitting Data
for order =1:9
    order
    [coeff, s] = polyfit(xx, yy, order);
    figure (order)
    plot(xx,yy, 'r.', 'MarkerSize', 24); hold on; 
    plot (xx, polyval(coeff,xx), 'k-', 'LineWidth',3);
    title (order)
    hold off; 
end 

%checking fitting
figure(19)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2= fit_metric.rsquare
    
    hold on 
    plot(order, r2, 'r.')
    
    title("r2");
    
end 
hold off 

figure(20)
for order =1:9
    order
    
    fitmodel = strcat('poly',num2str(order));
    [fit_out, fit_metric] = fit(xx, yy, fitmodel)
    r2adj = fit_metric.adjrsquare
    hold on 
    plot(order, r2adj, 'b.')
    title("r2adj")
    
end
hold off 

% we more data value, my analysis manages to find that the right order for
% the model. with only 10 data points, it seemed order 2 was already good
% fit, with more data, I can see that it is actaually order of 3 with
% better fit. 

% Problem 2. Basic statistics. 
% Part 1. Consider two different distributions - Gaussian numbers with a mean of
% 0 and variance 1 and Gaussian numbers with a mean of 1 and variance 1.
% (1) Make a plot of the average p-value for the t-test comparing N random
% numbers chosen from each of these two distributions as a function of N.
for N = 1:1000
xx = randn(N,1);
yy = 1+ sqrt(1) *randn(N,1);
[is_sig, pval] = ttest2(xx,yy);
hold on
plot(N, pval, 'r.');
end 
hold off 
xlabel("N"); ylabel("p-value"); 

% Part 2. Now keep the first distribution the same, but vary the mean of
% the second distribution between 0 and 10 with the same variance and
% repeat part one. Make a plot of all of these different curves on the same
% set of axes. What is special about the case where the mean of the second
% distribution is 0? 
for N = 1:100
    for a = 1:11
    xx = randn(N,1);
    yy = (a-1)+ sqrt(1) *randn(N,1);
    [is_sig(N,a), pval(N,a)] = ttest2(xx,yy);
   
    end
end 
plot(pval, 'o');
xlabel("N"); ylabel("p-value"); 
leg = 0:10;
leg = strtrim(cellstr(num2str(leg'))')
legend(leg); 
% when second distribution is 0, p-value is distributed randomly. 

% Part 3. Now keep the means of the two distributions at 0 and 1 as in part
% 1, but vary the variance of both distributions simultaneiously between 0.1 and 10 and plot the 
% p-values vs the number of numbers drawn as before. Comment on your results.  
for N = 1:100
    a = 0.1
    b = 1
    
    while a <=11
   
    xx = 0+ sqrt(a) *randn(N,1);
    yy = 1+ sqrt(a) *randn(N,1);
    [is_sig(N,b), pval(N,b)] = ttest2(xx,yy);
   b = b+1;
   a = a+1;
    end
end
%it's confusing when I plot everything on one plot -> I will make seperate
%plots
for i = 1:10
figure(i) 
plot(pval(:,i), 'o');
xlabel("N"); ylabel("p-value"); 
title(i)
end 

% one variance increases p-value is still relatively big even for greater
% N. 