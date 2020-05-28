from matplotlib.pyplot import figure, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
h = data[0]
e = data[1]

figure()
plot(h, e)
xlabel(r'mesh_size$h$')
ylabel(r'errors $e$')
savefig(output_file)

print('Generated ' + output_file)
