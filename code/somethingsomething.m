%%
data1_1 = cleandatami2(b25t450);
data1_2 = cleandatami2(b25t500);
data1_3 = cleandatami2(b25t450);

data2_1 = cleandatami2(b50t450);
data2_2 = cleandatami2(b50t500);
data2_3 = cleandatami2(b50t550);

data3_1 = cleandatami2(b75t450);
data3_2 = cleandatami2(b75t500);
data3_3 = cleandatami2(b75t550);


mines = [min(data1_1(:,3)) min(data1_2(:,3)) min(data1_3(:,3)) min(data2_1(:,3)) min(data2_2(:,3)) min(data2_3(:,3)) min(data3_1(:,3)) min(data3_2(:,3)) min(data3_3(:,3))]
mines = mines(end:-1:1)
[minmin indmin] = min(mines)

minconv = [data1_1(end,4) data1_2(end,4) data1_3(end,4) data2_1(end,4) data2_2(end,4) data2_3(end,4) data3_1(end,4) data3_2(end,4) data3_3(end,4)]
minconv = minconv(end:-1:1)
[minminconv indmonconv] = min(minconv)


