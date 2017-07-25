#!/usr/bin/octave --silent

dir = './post';

data = load(fullfile(dir, 'data.txt'));

h = figure;
hold off;
plot(data(:,1),data(:,2),"linewidth", 3, 'r');
hold on;
plot(data(:,1),data(:,3),"linewidth", 3, 'b');
hold on;
xlabel("time step","FontName","DejaVuSansMono","fontsize", 20);
ylabel("tracer input and output","FontName","DejaVuSansMono","fontsize", 20);
print(h, '-dsvg', fullfile(dir, "figure1.svg"));
refresh;
