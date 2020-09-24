import os
import numpy as np

#Set path names
path = os.getcwd()
data_path = os.path.join(os.path.split(path)[0],'files')
input_path = data_path+'\\example_reaction_list.txt'

f = open(input_path, "r")
f0 = f.read()
f1 = f0.split('\n')

rv = np.load( os.path.join(data_path,"rxn_vector.npy")).astype('str')
f2 = np.zeros(len(rv))
c = 0
for i in rv:
    if i+'0' in f1:
        f2[c] = 1
    c = c+1
test = []
for i in f1:
    if i[:-1] not in rv:
        test.append(i)

print(len(test))

np.save(os.path.join(data_path,'example_binary.npy'), f2)
