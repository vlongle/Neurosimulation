import numpy as np
import sys
import pandas as pd

def parse_data(filename, opt):
    with open(filename, 'r') as f:
       lines = f.readlines()

    lines = lines[3:] # eliminate 'Data generated at', '\n', '\n'

    endpoints = opt + 1
    cutoff = int(len(lines)/opt)

    begin_endpt = 0
    Xs = []


    # e.g. Xs = [data_1000, data_3000, data_5000]
    # e.g Xs = [data(1000, 1500), data(1500, 1000)]
    for i in np.linspace(0, len(lines), endpoints):
        if i == 0:
            continue

        end_endpt = int(i)

        line = [l.strip().split(' ') for l in lines[begin_endpt: end_endpt]]
        line = np.array(line, dtype=np.float64)
        #print('line', line.shape, line)
        Xs.append(line)
        begin_endpt = end_endpt


    # for subsequent data list in Xs. Only keep the last two columns (8, 9)
    # corresponding to fE, fI
    for i, X in enumerate(Xs[1:]):
        #print('X', X.shape)
        Xs[i+1] = X[:, [8,9]]

    # delete info about lambdaI, lambdaE, firingSpikesE, firingSpikeI
    Xs[0] = np.delete(Xs[0], range(4, 8), axis=1)

    return np.hstack(Xs)





if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Must pass in filename and int')
        exit(-1)

    filename = str(sys.argv[1])
    opt = int(sys.argv[2])

    if opt != 3 and opt != 2:
        print('Valid opt: 2 or 3')
        exit(-1)

    ret = parse_data(filename, opt)
    columns_name = ['SEE', 'SIE', 'SEI', 'SII'] + ['firing ' + str(l) for l\
            in range(0,len(ret[0])-4 )]
    df = pd.DataFrame(ret, columns=columns_name)
    print(df)

    name = filename.split('.')[0]
    np.savetxt(name+'_parsed.txt', ret)
