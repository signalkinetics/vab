function v = write_complex_binary(data, filename)

data_float = zeros(1, 2*length(data));
data_float(1, 1:2:end) = real(data);
data_float(1, 2:2:end) = imag(data);

v = write_float_binary(data_float, filename);