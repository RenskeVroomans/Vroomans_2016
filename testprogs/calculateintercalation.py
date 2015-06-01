import math

summ=0
p=0.05 #max prob of intercalation
d=0.075 #decay of morphogen
c=0 #nr of cells present at a certain time

for t in range(0, 120):
  c+=1
  inter=-p/d * math.exp(-d*c)+p/d
  summ+=inter
  c+=inter 
print "total sum is "+ str(summ) +"\ntotal nr of cells is "+str(c)
