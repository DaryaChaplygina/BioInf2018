f = open('mer_counts.txt', 'r')

int_ = ''
for line in f:
    if line.startswith('>'):
        int_ = line[1:]
    else:
        print(line.replace('\n', ''), int(int_))
