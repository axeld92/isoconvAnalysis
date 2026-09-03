function output = fs_mixture2(parameters,temp)

%calcula modelo de mezcla para tres componentes de Fraser suzuki



parameters_1 = parameters(1:4);
parameters_2 = parameters(5:8);
parameters_3 = parameters(9:12);



output = fs_function(temp, parameters_1) + fs_function(temp, parameters_2) + fs_function(temp, parameters_3);