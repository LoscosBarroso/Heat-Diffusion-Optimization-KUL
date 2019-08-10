M = [308.6 317.46 320.85 324.37 324.37 320.85 317.46 308.6 
303.13 314.12 319.12 322.62 322.62 319.12 314.12 303.13 
298.3 310.75 314.13 318 318 314.13 310.75 298.3 
293 303.17 308.17 313.91 313.91 308.17 303.17 293 
293 303.17 308.17 313.91 313.91 308.17 303.17 293 
298.3 310.75 314.13 318 318 314.13 310.75 298.3 
303.13 314.12 319.12 322.62 322.62 319.12 314.12 303.13 
308.6 317.46 320.85 324.37 324.37 320.85 317.46 308.6 
];
mean2(M)
figure;
heatmap(M);
v = [0 1 0.69362 0 0.69362 1 0 0 
0.56178 0.87251 0 0 0 0 0.87251 0.56178 
0 0 0 0 0 0 0 0 
0.85209 1 0.41072 0.42591 0.42591 0.41072 1 0.85209 
0 1 0.65955 0 0.65955 1 0 0 
1 1 0 0 0 0 1 1 
1 0.58857 0.38528 0 0.38528 0.58857 1 0 
0 0 0 0 0 0 0 0 
1 0.58857 0.38528 0 0.38528 0.58857 1 0 
1 1 0 0 0 0 1 1 
0 1 0.65955 0 0.65955 1 0 0 
0.85209 1 0.41072 0.42591 0.42591 0.41072 1 0.85209 
0 0 0 0 0 0 0 0 
0.56178 0.87251 0 0 0 0 0.87251 0.56178 
0 1 0.69362 0 0.69362 1 0 0 
];
figure;
set(gca,'visible','off')
xlim([0 1]);
ylim([0 1]);
plotVT(v);

saveas(gcf,'temps4.png');