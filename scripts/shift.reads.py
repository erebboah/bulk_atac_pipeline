import sys

infile = open(sys.argv[1],'w')

for line in sys.stdin:
	if line[0] == "@":
		infile.write(line.strip())
	else:
		f = line.strip().split('\t')
		if f[2] != '*':
			if f[1] == "99" or f[1] == "147":
				infile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                  (f[0], f[1], f[2], str(int(f[3]) + 4),f[4], f[5],f[6],f[7],f[8],f[9], f[10],f[11],f[12],f[13]))
			elif f[1] == "83" or f[1] == "163":
				infile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                  (f[0], f[1], f[2], str(int(f[3]) - 5),f[4], f[5],f[6],f[7],f[8],f[9], f[10],f[11],f[12],f[13]))
infile.close()
