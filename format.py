import sys
file1 = open("thyroid.txt")
i=0
count = 0
for line in file1:
    if i==-1:
        i+=1
        continue
    line = line.replace("\n","")
    line = line.split(" ")
    #if line[0] == "0" or line[0] == "3":
    for i in range(len(line)):
        if line[i]== '':
            continue
        elif i == len(line)-1:
    	       sys.stdout.write(line[i]+"\n")
        else:
    	       sys.stdout.write(line[i]+" ")
               count+=1
