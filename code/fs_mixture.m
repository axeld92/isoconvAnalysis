function output = fs_mixture(temp,parameters_1,parameters_2,parameters_3)

%calcula modelo de mezcla para tres componentes de Fraser suzuki

height_1 = parameters_1(1);
skew_1  = parameters_1(2);
position_1 = parameters_1(3);
width_1 = parameters_1(4);
height_2 = parameters_2(1);
skew_2 = parameters_2(2);
position_2 = parameters_2(3);
width_2 = parameters_2(4);
height_3 = parameters_3(1);
skew_3 = parameters_3(2);
position_3 = parameters_3(3);
width_3 = parameters_3(4);


output = fs_function(temp, parameters_1) + fs_function(temp, parameters_2) + fs_function(temp, parameters_3);