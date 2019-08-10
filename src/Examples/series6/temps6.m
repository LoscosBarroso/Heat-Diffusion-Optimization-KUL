M = [306.12 307.76 310.21 312.05 313.4 315.07 315.07 313.4 312.05 310.21 307.76 306.12 
303.81 306.23 310.17 311.51 313.21 314.36 314.36 313.21 311.51 310.17 306.23 303.81 
301.87 304.39 307.68 309.17 310.99 312.21 312.21 310.99 309.17 307.68 304.39 301.87 
299.99 302.23 305.85 307.36 308.86 310.45 310.45 308.86 307.36 305.85 302.23 299.99 
293 298.39 300.46 302.6 304.91 307.65 307.65 304.91 302.6 300.46 298.39 293 
293 295.42 297.79 300.18 302.71 305.68 305.68 302.71 300.18 297.79 295.42 293 
293 295.42 297.79 300.18 302.71 305.68 305.68 302.71 300.18 297.79 295.42 293 
293 298.39 300.46 302.6 304.91 307.65 307.65 304.91 302.6 300.46 298.39 293 
299.99 302.23 305.85 307.36 308.86 310.45 310.45 308.86 307.36 305.85 302.23 299.99 
301.87 304.39 307.68 309.17 310.99 312.21 312.21 310.99 309.17 307.68 304.39 301.87 
303.81 306.23 310.17 311.51 313.21 314.36 314.36 313.21 311.51 310.17 306.23 303.81 
306.12 307.76 310.21 312.05 313.4 315.07 315.07 313.4 312.05 310.21 307.76 306.12 
];
mean2(M)
figure;
heatmap(M);
v = [0 1 1 0.95186 0.60442 0 0.60442 0.95186 1 1 0 0 
0.51525 1 0 0 0 0 0 0 0 0 1 0.51525 
0 0 0 0 0 0 0 0 0 0 0 0 
0.80284 1 0.34292 0.36027 0.36794 0.37444 0.37444 0.36794 0.36027 0.34292 1 0.80284 
0 0 0.90541 0.57697 0 0 0 0.57697 0.90541 0 0 0 
1 1 1 0 0 0.58752 0.58752 0 0 1 1 1 
1 1 1 0.89979 0.75452 0 0.75452 0.89979 1 1 1 0 
1 1 0 0 0 0 0 0 0 0 1 1 
1 0.77737 0.66101 0.51799 0.33243 0 0.33243 0.51799 0.66101 0.77737 1 0 
1 0 0 0 0 0 0 0 0 0 0 1 
0.80595 0.72803 0.62639 0.49592 0.32123 0 0.32123 0.49592 0.62639 0.72803 0.80595 0 
0 0 0 0 0 0 0 0 0 0 0 0 
0.80595 0.72803 0.62639 0.49592 0.32123 0 0.32123 0.49592 0.62639 0.72803 0.80595 0 
1 0 0 0 0 0 0 0 0 0 0 1 
1 0.77737 0.66101 0.51799 0.33243 0 0.33243 0.51799 0.66101 0.77737 1 0 
1 1 0 0 0 0 0 0 0 0 1 1 
1 1 1 0.89979 0.75452 0 0.75452 0.89979 1 1 1 0 
1 1 1 0 0 0.58752 0.58752 0 0 1 1 1 
0 0 0.90541 0.57697 0 0 0 0.57697 0.90541 0 0 0 
0.80284 1 0.34292 0.36027 0.36794 0.37444 0.37444 0.36794 0.36027 0.34292 1 0.80284 
0 0 0 0 0 0 0 0 0 0 0 0 
0.51525 1 0 0 0 0 0 0 0 0 1 0.51525 
0 1 1 0.95186 0.60442 0 0.60442 0.95186 1 1 0 0 
];
figure;
set(gca,'visible','off')
xlim([0 1]);
ylim([0 1]);
plotVT(v);
saveas(gcf,'temps6.png');
