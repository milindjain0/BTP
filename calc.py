f = open("test.txt","r")
a = 5
i = 0
means = []
inputs = []
for line in f:
	line = line.replace("\n","")
	line = line.split()
	tmp = []
	for l in line:
		tmp.append(int(l))
	inputs.append(tmp)


for 
print inputs