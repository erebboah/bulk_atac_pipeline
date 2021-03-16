import sys

infile = open(sys.argv[1])
outfile= open(sys.argv[2], 'w')

counter = 1

for line in infile:
	name = "Peak" + str(counter)
	line = line.strip()
	outfile.write(line + '\t' + name + '\t' + "0" + '\t' + "+" + "\n")
	counter += 1
infile.close()
outfile.close()
