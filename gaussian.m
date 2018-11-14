clc;
clear;
close all;
%%
a = [4,5,6,9,3,2];
m = mean(a);
sd = std(a);
%--------- Gaussian PDF
x = -50:0.0001:50;
y = round(exp(-(x-m).^2/(2*sd^2)./(sd*sqrt(2*pi))), 4);

%%
plot(x,y)
xlabel('Frequency', 'FontSize',20)
ylabel('Y: PDF', 'FontSize',20)
title('Gaussian Distribution', 'FontSize',20)

hold on
stem( m, 1)

txt = '\leftarrow \mu = 4.83';
text(m, 0.5 ,txt, 'FontSize',20)

hold on
% I = find(y == 0.341);
% stem( x(I), y(I))
% txt = '\leftarrow \sigma';
% text( x(I), y(I)./2 ,txt, 'FontSize',20)

z = abs(x - (m-sd));
[z, I] = min(z);
stem( [m-sd, m+sd], [y(I), y(I)])
txt = '\mu - \sigma \rightarrow';
text(m-sd-9, y(I)./5, txt, 'FontSize',20)
txt = '\leftarrow \mu + \sigma';
text(m+sd, y(I)./5 ,txt, 'FontSize',20)



