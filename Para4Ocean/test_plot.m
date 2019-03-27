function test_plot(Ftitle, varargin)

figure, 
for i=3: 2:length(varargin)
    plot(varargin{1}, varargin{i},'LineWidth',2);
    hold on;
end
xlabel('Centerline distance (m)');ylabel('m');
legend(varargin(2:2:length(varargin)));
set(gca,'FontSize',15)
title(Ftitle);
end