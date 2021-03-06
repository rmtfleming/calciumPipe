function [fitresult, gof] = exponentialFit(tfall, xfall)
%CREATEFIT(TFALL,XFALL)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tfall
%      Y Output: xfall
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 07-Feb-2013 09:52:12


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( tfall, xfall );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf];
opts.StartPoint = [1 0.5 -1];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'real data', 'exponential fit' );
% Label axes
% xlabel( 'tfall' );
% ylabel( 'xfall' );
% hold on
% grid on
% print('-dtiff','exp_fit');

