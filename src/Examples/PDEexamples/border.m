M = [359.13 367.93 376.28 383.95 390.79 396.65 401.45 405.11 407.58 408.82 408.82 407.58 405.11 401.45 396.65 390.79 383.95 376.28 367.93 359.13 
350.14 410.25 454.86 489.68 517.28 539.04 555.77 567.95 575.92 579.85 579.85 575.92 567.95 555.77 539.04 517.28 489.68 454.86 410.25 350.14 
340.59 437.2 512.38 571.75 618.76 655.61 683.75 704.14 717.42 723.96 723.96 717.42 704.14 683.75 655.61 618.76 571.75 512.38 437.2 340.59 
330.26 454.71 554.86 635.31 699.53 750.03 788.63 816.59 834.78 843.75 843.75 834.78 816.59 788.63 750.03 699.53 635.31 554.86 454.71 330.26 
318.97 465.65 586.16 684.24 763.14 825.49 873.27 907.93 930.51 941.64 941.64 930.51 907.93 873.27 825.49 763.14 684.24 586.16 465.65 318.97 
306.59 471.91 609.04 721.47 812.46 884.66 940.16 980.51 1006.8 1019.8 1019.8 1006.8 980.51 940.16 884.66 812.46 721.47 609.04 471.91 306.59 
293 475.51 625.74 749.27 849.69 929.68 991.34 1036.3 1065.6 1080.1 1080.1 1065.6 1036.3 991.34 929.68 849.69 749.27 625.74 475.51 293 
293 480.51 638.29 769.33 876.48 962.17 1028.4 1076.7 1108.3 1123.9 1123.9 1108.3 1076.7 1028.4 962.17 876.48 769.33 638.29 480.51 293 
293 484.37 646.73 782.41 893.85 983.25 1052.5 1103.1 1136.2 1152.6 1152.6 1136.2 1103.1 1052.5 983.25 893.85 782.41 646.73 484.37 293 
293 486.4 650.96 788.88 902.41 993.63 1064.4 1116.1 1149.9 1166.7 1166.7 1149.9 1116.1 1064.4 993.63 902.41 788.88 650.96 486.4 293 
293 486.4 650.96 788.88 902.41 993.63 1064.4 1116.1 1149.9 1166.7 1166.7 1149.9 1116.1 1064.4 993.63 902.41 788.88 650.96 486.4 293 
293 484.37 646.73 782.41 893.85 983.25 1052.5 1103.1 1136.2 1152.6 1152.6 1136.2 1103.1 1052.5 983.25 893.85 782.41 646.73 484.37 293 
293 480.51 638.29 769.33 876.48 962.17 1028.4 1076.7 1108.3 1123.9 1123.9 1108.3 1076.7 1028.4 962.17 876.48 769.33 638.29 480.51 293 
293 475.51 625.74 749.27 849.69 929.68 991.34 1036.3 1065.6 1080.1 1080.1 1065.6 1036.3 991.34 929.68 849.69 749.27 625.74 475.51 293 
306.59 471.91 609.04 721.47 812.46 884.66 940.16 980.51 1006.8 1019.8 1019.8 1006.8 980.51 940.16 884.66 812.46 721.47 609.04 471.91 306.59 
318.97 465.65 586.16 684.24 763.14 825.49 873.27 907.93 930.51 941.64 941.64 930.51 907.93 873.27 825.49 763.14 684.24 586.16 465.65 318.97 
330.26 454.71 554.86 635.31 699.53 750.03 788.63 816.59 834.78 843.75 843.75 834.78 816.59 788.63 750.03 699.53 635.31 554.86 454.71 330.26 
340.59 437.2 512.38 571.75 618.76 655.61 683.75 704.14 717.42 723.96 723.96 717.42 704.14 683.75 655.61 618.76 571.75 512.38 437.2 340.59 
350.14 410.25 454.86 489.68 517.28 539.04 555.77 567.95 575.92 579.85 579.85 575.92 567.95 555.77 539.04 517.28 489.68 454.86 410.25 350.14 
359.13 367.93 376.28 383.95 390.79 396.65 401.45 405.11 407.58 408.82 408.82 407.58 405.11 401.45 396.65 390.79 383.95 376.28 367.93 359.13 
];
mean2(M)
figure;
heatmap(M);
v = [1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 0 
];
figure;
plotVT(v);
