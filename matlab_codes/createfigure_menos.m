function createfigure1(X1, YMatrix1)
%CREATEFIGURE1(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 08-Nov-2016 12:24:42

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',13);
%% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[60 480]);
%% Uncomment the following line to preserve the Y-limits of the axes
 %ylim(axes1,[0 6]);
grid(axes1,'on');
hold(axes1,'all');

% Create xlabel
xlabel('T (number of periods)','FontSize',16.5);

% Create ylabel
ylabel('Expected out-of-sample performance (x100)','FontSize',16.5);

% Create title
title('Expected out-of-sample performance','FontWeight','bold',...
    'FontSize',18);

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'Marker','.','Color',[1 0 0],'DisplayName','Optimal MV');
set(plot1(2),'Marker','.','Color',[0 0 1],...
    'DisplayName','Optimal only risky assets');
set(plot1(3),'Marker','o','LineStyle','-.','Color',[1 0 0],...
    'DisplayName','Optimal plug-in MV');
set(plot1(4),'Marker','o','Color',[0 0 1],...
    'DisplayName','Optimal plug-in only risky assets');
set(plot1(5),'LineWidth',2.5,'LineStyle','-.',...
    'Color',[0 0.498039215803146 0],...
    'DisplayName','Global minimum variance');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.618489583333332 0.172890733056708 0.255208333333333 0.30782802209616],...
    'FontSize',15);

