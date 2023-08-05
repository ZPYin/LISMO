clc;
close all;

%% Parameter Definition
sigma1 = 1;
sigma2 = 4;
length2 = 4;
x = -30:0.01:30;

func1 = @(x) exp(- (x - 0).^2 / (2 * sigma1^2));
func2 = @(x) (abs(x) - sigma2 / 2) .* (length2 + sigma2 / 2 - abs(x)) > 0;
y = conv(func1(x), func2(x), 'same') / sum(func2(x)) * 2;

figure; 
hold on;
plot(x, func1(x));
plot(x, func2(x));
plot(x, y);

%%
x1 = -10:0.03:10;
y1 = -10:0.03:10;
[X1, Y1] = meshgrid(x1, y1);
z1 = NaN(size(X1));
for iX = 1:length(x1)
    for iY = 1:length(y1)
        [~, iPos] = min(abs(x - sqrt(x1(iX).^2 + y1(iY).^2)));
        z1(iX, iY) = y(iPos);
    end
end

%% Display
figure;
p1 = pcolor(X1, Y1, z1);
p1.EdgeColor = 'none';
colormap('gray');
pbaspect([1, 1, 1]);
export_fig(gcf, 'laser_image.png', '-r600');
