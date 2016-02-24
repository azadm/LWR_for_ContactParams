function [time, matrix32, labels32] = ...
    convert_data_from_53_to_32 ( matrix53, labels53 )

time = matrix53(:,1);

order_vector = [42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 41, 40, ...
    39, 7, 8, 9, 10, 11, 12, 13, 1, 2, 3, 23, 24, 25, 26, 27, 28, 29];

% plus 1 because of the timestamp at the first column
matrix32 = matrix53(:,order_vector+1)';

labels32 = labels53(order_vector);

