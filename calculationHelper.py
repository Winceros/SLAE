fX = open("x_progonka.txt");
fMatrix = open("matrix.txt");

x = [float(lineX) for lineX in fX.readlines()]

i = 0
fMatrix.readline()
for row in fMatrix.readlines():
	res = 0
	j = 0
	[A,B] = row.split('=')		
	for aj in A.split():		
		print aj," ",x[i+j]
		res += float(aj)*x[i+j]
		j+=1		
	print "calc = ",res,"; truly_ans = ",B	
	i= i+ 1

raw_input()	