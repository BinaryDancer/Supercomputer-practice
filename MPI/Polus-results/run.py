with open("run.sh", "w") as f:
	for data in [(40, 'mini'), (120, 'small'), (400, 'medium'), (2000, 'large'), (4000, 'extra')]:
		i = 1
		while (i < data[0] and i < 65):
			for j in range(3):
				f.write("mpisubmit.pl -p {proc} {dataset}.o --stdout {dataset}{proc}_{it}\n".format(proc=i, dataset=data[1], it=j))
			i *= 2

