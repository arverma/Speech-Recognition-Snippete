clc;
close;
a = [4,5,6,9,3,2]; % Ludo | 1=4, 2=5, 3=6, 4=9, 5=3, 6=2 
b = [12,17]; %Coin flip | Head = 12, Tail = 17
c = [19,20,25,20,19,27,32,26,25,23]; % Cricket | Run per over(10) 

m(1) = mean(a);
m(2) = mean(b);
m(3) = mean(c);

sd(1) = std(a);
sd(2) = std(b);
sd(3) = std(c);

%--------- Gaussian PDF
x = -40:0.3:80;
y = zeros(3,length(x));
for i = 1:3
    y(i,:) = exp(-(x-m(i)).^2/(2*sd(i)^2)./(sd(i)*sqrt(2*pi)));
end


stem(x, y(1,:),'r')
xlabel('Frequency', 'FontSize',20)
ylabel('Y: PDF', 'FontSize',20)
title('Multi-Gaussian Distribution', 'FontSize',15)

hold on
stem(x, y(2,:), 'b')
hold on
stem(x, y(3,:),'g')

hold on
stem( m, [1,1,1],'b')
txt = {'\leftarrow \mu1','\leftarrow \mu2','\leftarrow \mu3' };
text(m, [0.5,0.5,0.5] ,txt, 'FontSize',18)