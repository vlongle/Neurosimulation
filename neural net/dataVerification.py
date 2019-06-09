import numpy as np


# also use lemma 4.3 to generate data
NE = 300
NI = 100

PEE = 0.15
PIE = 0.5
PEI = 0.5
PII = 0.4


M = 100


def CEE(SEE):
    return NE*PEE*SEE

def CII(SII):
    return NI*PII*SII

def CEI(SEI):
    return NI*PEI*SEI

def CIE(SIE):
    return NE*PEI*SIE

def calculate_Erate(SEE, SIE, SEI, SII, kickE, kickI):
    num = kickE*(M + CII(SII)) - kickI*CEI(SEI)

    print('num E_rate', num)

    denom = (M-CEE(SEE))*(M+CII(SII)) + (CEI(SEI)*CIE(SIE))

    print('denom E_rate', denom)


    return num/denom

def calculate_Irate(SEE, SIE, SEI, SII, kickE, kickI):
    num = kickI*(M + CEE(SEE)) - kickE*CIE(SIE)

    denom = (M-CEE(SEE))*(M+CII(SII)) + (CEI(SEI)*CIE(SIE))


    return num/denom


train_file = 'training_data.txt'


# x-train: 6 rates (E/I) for kickI, kickE 1000, 3000, 5000
x_train =  np.loadtxt(train_file, usecols=range(0,6), dtype=np.float32)


# SEE, SIE, SEI, SII
y_train = np.loadtxt(train_file, usecols=range(6,10), dtype=np.float32)

print('rates', x_train[0])
print('SEE, SIE, SEI, SII', y_train[0])


x_compare = []
for y in [y_train[0]]:
    x = []
    SEE, SIE, SEI, SII = y
    for i in range(3):
        kickE = 1000 + i*2000
        kickI = kickE
        x.append(calculate_Erate(SEE, SIE, SEI, SII, kickE, kickI))
        x.append(calculate_Irate(SEE, SIE, SEI, SII, kickE, kickI))

    x_compare.append(x)


print(x_compare[0])