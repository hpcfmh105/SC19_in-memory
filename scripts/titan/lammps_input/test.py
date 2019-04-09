import os.path
path = "./"
dirs = os.listdir(path)
for file in dirs:
    if 'in.lj.32v16_' in file:
        f = open(file, 'r')
        w_str = ''
        for line in f:
            index = file
            print(index)
            index = index.strip('in.lj.32v16')
            index = index.strip('_')
            index = index.strip('.txt')
            print(index)
            if 'z index' in line:
                line = line.replace('number', index) 
                w_str += line
                print(line)
            else:
                w_str += line
        wopen=open(file,'w')
        wopen.write(w_str)
        f.close()
        wopen.close()
