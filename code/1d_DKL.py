import math
import numpy
import scipy
import itertools
from scipy.special import digamma
import matplotlib.pyplot as plt

class Point:
    def __init__(self, x):
        self._x = x


class Pair:
    def __init__(self, num, i):
        self._number = num
        self._index = i


def in_point(fname):
    infile = open(fname, "r")
    lst = infile.readlines()
    point_l = []
    infile.close()
    for line in lst:
        l = line.split()
        if ',' in l[0]:
            x = l[0][0] + l[0][2:]
            # y = float(l[1])
        # elif ',' in l[1]:
        #     x = float(l[0])
        #     y = l[1][0] + l[1][2:]
        else:
            x = float(l[0])
            # print(x)
            # y = float(l[1])
        point_l.append(Point(float(x)))
    return point_l


# square
# point_l = (Listof Point)

def prob(fname, h):
    point_l = in_point(fname)
    L1 = []
    k = 1
    n_last = 0
    N = len(point_l)
    list_x = []
    # list_y = []
    for point in point_l:
        list_x.append(point._x)
        # list_y.append(point._y)
    max_x = max(list_x)
    r=max_x/h
    # print(max_x)
    # max_y = max(list_y)
    # field = max(max_x, max_y)
    while (k - 1) * r <= max_x:
        n = 0
        for point in point_l:
            if (-k * r) <= point._x <= k * r :
                n += 1
        p = (n - n_last) / N
        L1.append(p)
        k += 1
        n_last = n
    return L1

def window_prob(fname,h,i,lg):
    point_l = in_point(fname)
    l = point_l[i:lg+i]
    L = []
    k = 1
    n_last = 0
    list_x = []
    # list_y = []
    for point in l:
        list_x.append(point._x)
        # list_y.append(point._y)
    #print(list_x)
    #print(list_y)
    max_x = max(list_x)
    r=max_x/h
    # max_y = max(list_y)
    # field = max(max_x,max_y)
    while (k-1)*r <= max_x:
        n = 0
        for point in l:
            if (-k*r) <= point._x <= k*r:
                n += 1
        p = (n - n_last) / lg
        L.append(p)
        k += 1
        n_last = n
    return L


def DKL(d1, d2, h):
    result = 0
    L1 = prob(d1, h)
    L2 = prob(d2, h)
    n = min(len(L1), len(L2))
    for i in range(n):
        if L1[i] != 0 and L2[i] != 0:
            result += L1[i] * (math.log(L2[i] / L1[i]))
    return result * (-1)


def DJS(d1, d2, h):
    result = 0
    M = []
    L1 = prob(d1, h)
    L2 = prob(d2, h)
    n = min(len(L1), len(L2))
    for i in range(n):
        M.append((L1[i] + L2[i]) / 2)
    if n == len(L1):
        for p in L2[n:]:
            M.append(p / 2)
    elif n == len(L2):
        for p in L1[n:]:
            M.append(p / 2)
    for i in range(n):
        if L1[i] != 0:
            result -= L1[i] * math.log(M[i] / L1[i])
        if L2[i] != 0:
            result -= L2[i] * math.log(M[i] / L2[i])
    return result / 2


def window_DJS(fname_1, fname_2, h, lg):
    x = []
    y = []
    if isinstance(fname_1, list):
        l = len(fname_1)
    else:
        l = max(len(in_point(fname_1)), len(in_point(fname_2)))
    for num in range(0, (l - lg), lg):
        result = 0
        M = []
        L1 = window_prob(fname_1, h, num, lg)
        L2 = window_prob(fname_2, h, num, lg)
        n = min(len(L1), len(L2))
        for i in range(n):
            M.append((L1[i] + L2[i]) / 2)
        if n == len(L1):
            for p in L2[n:]:
                M.append(p / 2)
        elif n == len(L2):
            for p in L1[n:]:
                M.append(p / 2)
        for i in range(n):
            if L1[i] != 0:
                result -= L1[i] * math.log(M[i] / L1[i])
            if L2[i] != 0:
                result -= L2[i] * math.log(M[i] / L2[i])
        result = result / 2
        x.append(num)
        y.append(result)
    # print(x)
    # print(y)
    y1 = []
    y2 = []
    for i in range(0, (l - lg), lg):
        y1.append(0)
        y2.append(1)
    # x = np.array(x)
    # y = np.array(y)
    plt.plot(x, y)
    plt.plot(x, y1, '--')
    plt.plot(x, y2, '--')
    plt.ylabel('window_DJS')
    plt.xlabel('window')
    plt.show()


def DJSa(d1, d2, h, a):
    assert a > -1 and a < 1, 'result is nan'
    result = 0
    M = []
    if isinstance(d1, list) and isinstance(d2, list):
        L1 = d1
        L2 = d2
    else:
        L1 = prob(d1, h)
        L2 = prob(d2, h)
    n = min(len(L1), len(L2))
    for i in range(n):
        M.append((L1[i] + L2[i]) / 2)
    if n == len(L1):
        for p in L2[n:]:
            M.append(p / 2)
    elif n == len(L2):
        for p in L1[n:]:
            M.append(p / 2)
    gam = math.gamma(a + 1)
    # print(gam)
    digam = digamma(1) - digamma(1 - a)
    # print(digam)
    for num in L1:
        if num != 0.0:
            result -= num * (num ** (-a)) * (math.log(num) + digam)
            # print(result,'@')
    for num in L2:
        if num != 0.0:
            result -= num * (num ** (-a)) * (math.log(num) + digam)
            # print(result,'#')
    result = result / 2
    # print(result)
    for num in M:
        if num != 0.0:
            result += num * ((num ** (-a))) * (math.log(num) + digam)
    return result / gam


def window_DJSa(fname_1, fname_2, h, a, lg):
    assert a > -1 and a < 1, 'result is nan'
    x = []
    y = []
    if isinstance(fname_1, list):
        l = len(fname_1)
    else:
        l = len(in_point(fname_1))
    for k in range(0, (l - lg), lg):
        result = 0
        M = []
        L1 = window_prob(fname_1, h, k, lg)
        L2 = window_prob(fname_2, h, k, lg)
        n = min(len(L1), len(L2))
        for i in range(n):
            M.append((L1[i] + L2[i]) / 2)
        if n == len(L1):
            for p in L2[n:]:
                M.append(p / 2)
        elif n == len(L2):
            for p in L1[n:]:
                M.append(p / 2)
        gam = math.gamma(a + 1)
        # print(gam)
        digam = digamma(1) - digamma(1 - a)
        # print(digam)
        for num in L1:
            if num != 0.0:
                result -= num * (num ** (-a)) * (math.log(num) + digam)
                # print(result,'@')
        for num in L2:
            if num != 0.0:
                result -= num * (num ** (-a)) * (math.log(num) + digam)
                # print(result,'#')
        result = result / 2
        # print(result)
        for num in M:
            if num != 0.0:
                result += num * ((num ** (-a))) * (math.log(num) + digam)
        x.append(k)
        y.append(result / gam)
    print(x)
    print(y)
    y1 = []
    y2 = []
    for i in range(0, (l - lg), lg):
        y1.append(0)
        y2.append(1)
    # print(y1)
    # print(y2)
    # x = np.array(x)
    # y = np.array(y)
    plt.plot(x, y)
    plt.plot(x, y1, '--')
    plt.plot(x, y2, '--')
    plt.ylabel('window_DJSa')
    plt.xlabel('window')
    plt.show()


def quicksort(pairs):
    if len(pairs) <= 1:
        return pairs
    less = []
    greater = []
    base = pairs.pop()
    base_num = base._number
    for x in pairs:
        if x._number < base_num:
            less.append(x)
        else:
            greater.append(x)
    return quicksort(less) + [base] + quicksort(greater)


def prob_pentropy(fname, n):
    assert n >= 2 and n <= 7, 'result is not accurate'
    i = 0
    order_l = []
    standard_l = {}
    point_l = in_point(fname)
    list_x = []
    # list_y = []
    for point in point_l:
        list_x.append(point._x)
        # list_y.append(point._y)
    lst = list_x
    while i + n <= len(lst):    # ori:i + n < len(lst)
        pair_l = []
        sub_lst = lst[i:i + n]
        k = 0
        for num in sub_lst:
            pair_l.append(Pair(num, str(k)))
            k += 1
        pair_l = quicksort(pair_l)
        order_str = ''
        for pair in pair_l:
            order_str += pair._index
        order_l.append(order_str)
        i += 1
    # print(order_l)
    permutations = list(itertools.permutations(range(n)))
    # print(permutations)
    for lst in permutations:
        std = ''
        for num in lst:
            std += str(num)
        standard_l[std] = 0
    # print('@',standard_l)
    for num in order_l:
        for std in standard_l:
            if num == std:
                standard_l[std] += 1
    # print('#',standard_l)
    for num in standard_l.keys():
        standard_l[num] = standard_l[num] / len(order_l)
    return standard_l


def window_prob_pentropy(fname, i, n, lg):  # data=0 or 1->x or y,lg->window length
    assert n >= 2 and n <= 7, 'result is not accurate'
    # i = 0
    point_l = in_point(fname)[i:lg + i]
    i = 0
    order_l = []
    standard_l = {}
    list_x = []
    # list_y = []
    for point in point_l:
        list_x.append(point._x)
        # list_y.append(point._y)
    lst = list_x
    while i + n <= len(lst):
        pair_l = []
        sub_lst = lst[i:i + n]
        k = 0
        for num in sub_lst:
            pair_l.append(Pair(num, str(k)))
            k += 1
        pair_l = quicksort(pair_l)
        order_str = ''
        for pair in pair_l:
            order_str += pair._index
        order_l.append(order_str)
        i += 1                    #ori: i +=1 for moving window
    # print(order_l)
    permutations = list(itertools.permutations(range(n)))
    for lst in permutations:
        std = ''
        for num in lst:
            std += str(num)
        standard_l[std] = 0
    # print('@',standard_l)
    for num in order_l:
        for std in standard_l:
            if num == std:
                standard_l[std] += 1
    # print('#',standard_l)
    for num in standard_l.keys():
        standard_l[num] = standard_l[num] / len(order_l)
    return standard_l


def DJS_PE(fname_1, fname_2, n):
    d1 = prob_pentropy(fname_1, n)
    d2 = prob_pentropy(fname_2, n)
    M = d1
    for num in d2.keys():
        M[num] = (M[num] + d2[num]) / 2
    # print(M)
    prob_d1 = list(filter(lambda x: x != 0, list(map(lambda k: d1[k], d1.keys()))))
    prob_d2 = list(filter(lambda x: x != 0, list(map(lambda k: d2[k], d2.keys()))))
    prob_m = list(filter(lambda x: x != 0, list(map(lambda k: M[k], M.keys()))))
    # print(prob_d1)
    pe_1 = -sum(prob_d1 * (numpy.log(prob_m)-numpy.log(prob_d1)))
    pe_2 = -sum(prob_d2 * (numpy.log(prob_m)-numpy.log(prob_d2)))
    # pe_m = -sum(prob_m * numpy.log(prob_m))
    result = 1 / 2 * (pe_1 + pe_2)
    return result

def DJSa_PE(d1,d2,a,n):
    assert a > -1 and a < 1, 'result is nan'
    result = 0
    L1 = prob_pentropy(d1,n)
    L2 = prob_pentropy(d2,n)
    M = L1
    gam = math.gamma(a+1)
    #print(gam)
    digam = digamma(1) - digamma(1 - a)
    #print(digam)
    for num in L2.keys():
        M[num] = (M[num] + L2[num]) / 2
    prob_d1 = list(filter(lambda x: x!=0, list(map(lambda k: L1[k],L1.keys()))))
    prob_d2 = list(filter(lambda x: x!=0, list(map(lambda k: L2[k],L2.keys()))))
    prob_m = list(filter(lambda x: x!=0, list(map(lambda k: M[k], M.keys()))))
    for num in prob_d1:
        if num != 0.0:
            result -= num * (num**(-a)) * (math.log(num) + digam)
            #print(result,'@')
    for num in prob_d2:
        if num != 0.0:
            result -= num * (num**(-a)) * (math.log(num) + digam)
            #print(result,'#')
    result = result / 2
    #print(result)
    for num in prob_m:
        if num != 0.0:
            result += num*((num**(-a))) * (math.log(num) + digam)
    return result / gam


def window_DJS_PE(fname_1, fname_2, n, lg):
    X = []
    Y = []
    l = max(len(in_point(fname_1)), len(in_point(fname_2)))
    for k in range(0, (l - lg), lg):
        result = 0
        M = []
        d1 = window_prob_pentropy(fname_1, k, n, lg)
        d2 = window_prob_pentropy(fname_2, k, n, lg)
        M = d1
        for num in d2.keys():
            M[num] = (M[num] + d2[num]) / 2
        # print(M)
        # print("--------------------------")
        # for i in range(0, len(M) - k): #not sure for moving window!
        prob_d1 = list(filter(lambda x: x != 0, list(map(lambda k: d1[k], d1.keys()))))
        prob_d2 = list(filter(lambda x: x != 0, list(map(lambda k: d2[k], d2.keys()))))
        prob_m = list(filter(lambda x: x != 0, list(map(lambda k: M[k], M.keys()))))
        pe_1 = -sum(prob_d1 * (numpy.log(prob_m) - numpy.log(prob_d1)))
        pe_2 = -sum(prob_d2 * (numpy.log(prob_m) - numpy.log(prob_d2)))
        result = 1 / 2 * (pe_1 + pe_2)
        # print(result)
        # print("--------------------------")
        X.append(k)
        Y.append(result)
    # print(X)
    # print(Y)
    y1 = []
    y2 = []
    # for i in range(0, (l - lg), lg):
    #     y1.append(0)
    #     y2.append(1)
    # x = np.array(x)
    # y = np.array(y)
    plt.plot(X, Y)
    # plt.plot(X, y1, '--')
    # plt.plot(X, y2, '--')
    plt.ylabel('window_DJS_PE')
    plt.xlabel('window')
    plt.show()

import math
import numpy as np
import scipy
import itertools
#from scipy.special import digamma
#import matplotlib.pyplot as plt

def a(n,d):
    x = d * math.gamma(n - d) / (math.gamma(1 - d) * math.gamma(n + 1))
    return x
def a1(n,d):
    #striling推导的计算gamma函数方法
    x = (1-(d+1)/(n+1))**(n+1/2)*math.exp(1+d)*d/((n-d)**(d+1))/math.gamma(1-d)
    return x
def afm(d1,W):
    # print '%d %d %d' %(n-d,1-d,n+1)
    N=5000;
    np.random.seed(1)
    mu1, sigma1, n_samples1 = 0, 0.1, N  # mean, standard deviation, size
    mu2, sigma2, n_samples2 = 0, 0.1, N
    w1 = np.random.normal(mu1, sigma1, n_samples1)
    w2 = np.random.normal(mu2, sigma2, n_samples2)
    # print(w1,w2)
    # return result
    d2=d1
    X=[0]*N
    Y=[0]*N
    # X[1], Y[1]=0, 0
    # W = 0.5;
    x = [0] * N
    y = [0] * N
    x[0]=w1[0]
    y[0]=w2[0]
    list_p=np.ones([N,2])
    # d1 = 0.3
    # d2 = 0.2
    # t = 0
    list_p[0][0]=w1[0]
    list_p[0][1]=w2[0]
    for t in range(1,N):
        n = 1
        while n >= 1 and n <= t:
            X[t] += x[t - n] * a1(n, d1)
            # print(a1(n, d1))
            Y[t] += y[t - n] * a1(n, d2)
            # print('%f---%f' % (X[t], Y[t]))
            n = n + 1
        x[t] = W * X[t] + (1 - W) * Y[t] + w1[t]
        y[t] = (1 - W) * X[t] + W * Y[t] + w2[t]
        list_p[t][0] = x[t]
        list_p[t][1] = y[t]
        # print('%f---%f' % (x[t], y[t]))
    # return x,y
    # plt.plot(x, '-')
    # plt.plot(y, '-')
    # plt.show()
    # print(X,Y)
    # print(1)
    # print(list_p)
    # np.savetxt('simulation', list_p)
    np.savetxt('simulation1', x)
    np.savetxt('simulation2', y)
    return list_p
