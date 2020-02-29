import math
import numpy
import scipy
import itertools
from scipy.special import digamma
import matplotlib.pyplot as plt
import xlwt


class Point:
    def __init__(self,x,y):
        self._x = x
        self._y = y

class Pair:
    def __init__(self,num,i):
        self._number = num
        self._index = i
        
# def in_point(fname): #for txt
#     infile = open(fname, "r")
#     lst = infile.readlines()
#     point_l = []
#     infile.close()
#     for line in lst:
#         l = line.split()
#         if ',' in l[0]:
#             x = l[0][0]+l[0][2:]
#             y = float(l[1])
#         elif ',' in l[1]:
#             x = float(l[0])
#             y = l[1][0]+l[1][2:]
#         else:
#             x = float(l[0])
#             y = float(l[1])
#         point_l.append(Point(float(x),float(y)))
#     return point_l

def in_point(fname):  #for csv
    nd = numpy.genfromtxt(fname, delimiter=',', skip_header=True)
    final_list = nd.tolist()
    point_l=[]
    for l in final_list:
        x = float(l[5])
        y = float(l[6])
        if not(math.isnan(x)) and not(math.isnan(y)) and y!=0.0:
            point_l.append(Point(float(x), float(y))) #Point
    return point_l

def in_point2(fname):  #for csv
    nd = numpy.genfromtxt(fname, delimiter=',', skip_header=True)
    final_list = nd.tolist()
    point_l=[]
    for l in final_list:
        # print(l)
        x = float(l[5])
        # y = float(l[6])
        point_l.append(float(x)) #Point
    return point_l

# square
# point_l = (Listof Point)

def prob1(fname,h):   #rectangle
    point_l = in_point(fname)
    L1 = []
    k = 1
    n_last = 0
    N = len(point_l)
    list_x = []
    list_y = []
    for point in point_l:
        list_x.append(point._x)
        list_y.append(point._y)
    max_x = max(list_x)
    max_y = max(list_y)
    # field = max(max_x/h,max_y/h)
    r1 = max_x / h
    r2 = max_y / h
    while k-1<= h:
        n = 0
        for point in point_l:
            if (-k * r1) <= point._x <= k * r1 and (-k * r2) <= point._y <= k * r2:
                n += 1
        p = (n - n_last) / N
        L1.append(p)
        k += 1
        n_last = n
    return L1

def prob_fail2(fname, h): #simplified version
    point_l = fname
    L1=[]
    for i in range(2):
        if point_l[0]!= point_l[1]:
            point_l[i]=point_l[i]-min(point_l)
        else:
            point_l=[0.5, 0.5]
    # print('@@@', point_l)
    for i in range(2):
        L1.append(point_l[i]/(sum(point_l)+1/1000000))
    return L1

def prob_fail(fname, h): #sosososimplified version
    point_l = fname
    L1=[]
    for i in range(2):
        L1.append(point_l[i]+0.5)
    return L1


def prob(d):   #square
    L1 = []
    for k in d:
        L1.append((k-min(d))/sum(d))
    return L1

def window_prob0(fname,r,i,lg):  #original version
    point_l = in_point(fname)
    l = point_l[i:lg+i]
    L = []
    k = 1
    n_last = 0
    list_x = []
    list_y = []
    for point in l:
        list_x.append(point._x)
        list_y.append(point._y)
    #print(list_x)
    #print(list_y)
    max_x = max(list_x)
    max_y = max(list_y)
    field = max(max_x,max_y)
    while (k-1)*r <= field:
        n = 0
        for point in l:
            if (-k*r) <= point._x <= k*r and (-k*r)<= point._y <= k*r:
                n += 1
        p = (n - n_last) / lg
        L.append(p)
        k += 1
        n_last = n
    return L

def window_prob(fname,h,i,lg): #square
    point_l = in_point(fname)
    l = point_l[i:lg+i]
    L = []
    k = 1
    n_last = 0
    list_x = []
    list_y = []
    for point in l:
        list_x.append(point._x)
        list_y.append(point._y)
    #print(list_x)
    #print(list_y)
    max_x = max(list_x)
    max_y = max(list_y)
    field = max(max_x,max_y)
    r = field / h
    while k-1<= h:
        n = 0
        for point in l:
            if (-k*r) <= point._x <= k*r and (-k*r)<= point._y <= k*r:
                n += 1
        p = (n - n_last) / lg
        L.append(p)
        k += 1
        n_last = n
    return L

def DKL(d1,d2,h):
    result = 0
    if isinstance(d1,list) and isinstance(d2,list):
        L1 = d1
        L2 = d2
    else:
        L1 = prob(d1,h)
        L2 = prob(d2,h)
    n = min(len(L1) ,len(L2))   
    for i in range(n):
        if L1[i] != 0 and L2[i] != 0:
            result -= L1[i]*(math.log(L2[i] / L1[i]))
    return result

    
def DJS(d1,d2):
    result = 0
    M = []
    L1 = prob([d1.tolist()[0],d2.tolist()[0]])
    L2 = prob([d1.tolist()[1],d2.tolist()[1]])
    n = min(len(L1) ,len(L2))   
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
            result -= L1[i]*math.log(M[i] / L1[i])
        if L2[i] != 0:
            result -= L2[i]*math.log(M[i] / L2[i])    
    return result / 2

def window_DJS(fname_1,fname_2,h,lg):
    x = []
    y = []
    l = max(len(in_point(fname_1)),len(in_point(fname_2)))
    for num in range(0,(l-lg), lg):
        result = 0
        M = []
        L1 = window_prob(fname_1,h,num,lg)
        L2 = window_prob(fname_2,h,num,lg)
        n = min(len(L1) ,len(L2))   
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
                result -= L1[i]*math.log(M[i] / L1[i])
            if L2[i] != 0:
                result -= L2[i]*math.log(M[i] / L2[i])    
        result = result / 2
        x.append(num)
        y.append(result)
    #print(x)
    #print(y)
    y1 = []
    y2 = []
    for i in range(0,(l-lg), lg):
        y1.append(0)
        y2.append(1)
    #x = np.array(x)
    #y = np.array(y)
    plt.plot(x,y)
    plt.plot(x,y1,'--')
    plt.plot(x,y2,'--')
    plt.ylabel('window_DJS')
    plt.xlabel('window')
    plt.show()

 
    
def DJSa(d1,d2,h,a):
    assert a > -1 and a < 1, 'result is nan'
    result = 0
    M = []
    L1 = prob(d1,h)
    L2 = prob(d2,h)
    n = min(len(L1) ,len(L2))   
    for i in range(n):
        M.append((L1[i] + L2[i]) / 2)
    if n == len(L1):
        for p in L2[n:]:
            M.append(p / 2)
    elif n == len(L2):
        for p in L1[n:]:
            M.append(p / 2)
    gam = math.gamma(a+1)
    #print(gam)
    digam = digamma(1) - digamma(1 - a)
    #print(digam)
    for num in L1:
        if num != 0.0:
            result -= num * (num**(-a)) * (math.log(num) + digam)
            #print(result,'@')
    for num in L2:
        if num != 0.0:
            result -= num * (num**(-a)) * (math.log(num) + digam)
            #print(result,'#')
    result = result / 2
    #print(result)
    for num in M:
        if num != 0.0:
            result += num*((num**(-a))) * (math.log(num) + digam)
    return result / gam  



def window_DJSa(fname_1, fname_2, h, a, lg):
    assert a > -1 and a < 1, 'result is nan'
    x = []
    y = []
    l = len(in_point(fname_1))
    for k in range(0,(l-lg), lg):
        result = 0
        M = []
        L1 = window_prob(fname_1,h,k,lg)
        L2 = window_prob(fname_2,h,k,lg)
        n = min(len(L1) ,len(L2))   
        for i in range(n):
            M.append((L1[i] + L2[i]) / 2)
        if n == len(L1):
            for p in L2[n:]:
                M.append(p / 2)
        elif n == len(L2):
            for p in L1[n:]:
                M.append(p / 2)
        gam = math.gamma(a+1)
        #print(gam)
        digam = digamma(1) - digamma(1 - a)
        #print(digam)
        for num in L1:
            if num != 0.0:
                result -= num * (num**(-a)) * (math.log(num) + digam)
                #print(result,'@')
        for num in L2:
            if num != 0.0:
                result -= num * (num**(-a)) * (math.log(num) + digam)
                #print(result,'#')
        result = result / 2
        #print(result)
        for num in M:
            if num != 0.0:
                result += num*((num**(-a))) * (math.log(num) + digam)
        x.append(k)
        y.append(result / gam) 
    # print(x)
    # print(y)
    y1 = []
    y2 = []
    for i in range(0,(l-lg), lg):
        y1.append(0)
        y2.append(1)
    #x = np.array(x)
    #y = np.array(y)
    plt.plot(x,y)
    plt.plot(x,y1,'--')
    plt.plot(x,y2,'--')
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

def prob_pentropy(fname,n,data):
    assert n >= 2 and n <= 7, 'result is not accurate'
    i = 0
    order_l = []
    standard_l = {}
    point_l = in_point(fname) 
    list_x = []
    list_y = []
    for point in point_l:
        list_x.append(point._x)
        list_y.append(point._y) 
    if data == 1:
        lst = list_x
    else:
        lst = list_y
    while i+n < len(lst):
        pair_l = []
        sub_lst = lst[i:i+n]
        k = 0
        for num in sub_lst:
            pair_l.append(Pair(num,str(k)))
            k += 1
        pair_l = quicksort(pair_l)
        order_str = ''
        for pair in pair_l:
            order_str += pair._index
        order_l.append(order_str)
        i += 1
    #print(order_l)   
    permutations = list(itertools.permutations(range(n)))
    for lst in permutations:
        std =''
        for num in lst:
            std += str(num)
        standard_l[std] = 0
    #print('@',standard_l)
    for num in order_l:
        for std in standard_l:
            if num == std:
                standard_l[std] += 1
    #print('#',standard_l)
    for num in standard_l.keys():
        standard_l[num] = standard_l[num] / len(order_l)
    return standard_l

def window_prob_pentropy(fname,i,n,data,lg): #data=0 or 1->x or y,lg->window length
    assert n >= 2 and n <= 7, 'result is not accurate'
    point_l = in_point(fname)[i:lg+i]
    i = 0
    order_l = []
    standard_l = {}
    list_x = []
    list_y = []
    for point in point_l:
        list_x.append(point._x)
        list_y.append(point._y) 
    if data == 1:
        lst = list_x
    else:
        lst = list_y
    while i+n < len(lst):
        pair_l = []
        sub_lst = lst[i:i+n]
        k = 0
        for num in sub_lst:
            pair_l.append(Pair(num,str(k)))
            k += 1
        pair_l = quicksort(pair_l)
        order_str = ''
        for pair in pair_l:
            order_str += pair._index
        order_l.append(order_str)
        i += 1
    #print(order_l)   
    permutations = list(itertools.permutations(range(n)))
    for lst in permutations:
        std =''
        for num in lst:
            std += str(num)
        standard_l[std] = 0
    #print('@',standard_l)
    for num in order_l:
        for std in standard_l:
            if num == std:
                standard_l[std] += 1
    #print('#',standard_l)
    for num in standard_l.keys():
        standard_l[num] = standard_l[num] / len(order_l)
    return standard_l
  
    
def DJS_PE(fname_1,fname_2,n,data):
    d1 = prob_pentropy(fname_1,n,data)
    d2 = prob_pentropy(fname_2,n,data)
    M = d1   
    for num in d2.keys():
        M[num] = (M[num] + d2[num]) / 2
    # print(M)
    prob_d1 = list(filter(lambda x: x!=0, list(map(lambda k: d1[k],d1.keys()))))
    prob_d2 = list(filter(lambda x: x!=0, list(map(lambda k: d2[k],d2.keys()))))
    prob_m = list(filter(lambda x: x!=0, list(map(lambda k: M[k], M.keys()))))
    pe_1 = -sum(prob_d1 * (numpy.log(prob_m) - numpy.log(prob_d1)))
    pe_2 = -sum(prob_d2 * (numpy.log(prob_m) - numpy.log(prob_d2)))
    # pe_m = -sum(prob_m * numpy.log(prob_m))
    result = 1 / 2 * (pe_1 + pe_2)
    return result 

def DJSa_PE(d1,d2,a,n,data):
    assert a > -1 and a < 1, 'result is nan'
    result = 0
    L1 = prob_pentropy(d1,n,data)
    L2 = prob_pentropy(d2,n,data)
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

def window_DJS_PE(fname_1,fname_2,n,data,lg):
    X = []
    Y = []
    if isinstance(fname_1,list):
        l = max(len(fname_1),len(fname_2))
    else:
        l = max(len(in_point(fname_1)),len(in_point(fname_2)))
    for k in range(0,(l-lg), lg):
        M = []    
        d1 = window_prob_pentropy(fname_1, k, n, data,lg)
        d2 = window_prob_pentropy(fname_2, k, n, data,lg)
        M = d1
        for num in d2.keys():
            M[num] = (M[num] + d2[num]) / 2
        # print(M)
        # for i in range(0, len(M) - k): #not sure for moving window!
        prob_d1 = list(filter(lambda x: x!=0, list(map(lambda k: d1[k],d1.keys()))))
        prob_d2 = list(filter(lambda x: x!=0, list(map(lambda k: d2[k],d2.keys()))))
        prob_m = list(filter(lambda x: x!=0, list(map(lambda k: M[k], M.keys()))))
        pe_1 = -sum(prob_d1 * (numpy.log(prob_m) - numpy.log(prob_d1)))
        pe_2 = -sum(prob_d2 * (numpy.log(prob_m) - numpy.log(prob_d2)))
        result = 1 / 2 * (pe_1 + pe_2)
        X.append(k)
        Y.append(result)
    # print(X)
    # print(Y)
    # y1 = []
    # y2 = []
    # for i in range(0,(l-lg),lg):
    #     y1.append(0)
    #     y2.append(1)
    #x = np.array(x)
    #y = np.array(y)
    plt.plot(X,Y)
    # plt.plot(X,y1,'--')
    # plt.plot(X,y2,'--')
    plt.ylabel('window_DJS_PE')
    plt.xlabel('window')
    plt.show()

class MakeMatrix_DJS:
    def inputinfor(self):
        while True:
            h=int(input("h "))
            n=int(input("n "))
            q=[]
            d=[[] for i in range(n)]
            for i in range(n):
                q.append(prob(input("fname"), h))
            # print('@', q)

            for i in range(n):
                for j in range(n):
                    n0 = min(len(q[i]), len(q[j]))
                    # print('@', n0)
                    result = 0
                    for k in range(n0):
                        # print(q[i][k],'---',q[j][k])
                        if q[i][k] != 0:
                            result -= q[i][k] * math.log(((q[i][k] + q[j][k]) / 2) / q[i][k])
                        if q[j][k] != 0:
                            result -= q[j][k] * math.log(((q[i][k] + q[j][k]) / 2) / q[j][k])
                    d[i].append(result / 2)
            # return d
            # print('@', d)
            #         print('@', result/2)
            #         print("--------------------")
            #     print("--------------------")
            # print("--------------------")
            def set_style(name, height, bold=False):
                style = xlwt.XFStyle()  # 初始化样式
                font = xlwt.Font()  # 为样式创建字体
                font.name = name
                font.bold = bold
                font.color_index = 4
                font.height = height
                style.font = font
                return style

            def write_excel():
                # 创建工作簿
                workbook = xlwt.Workbook(encoding='utf-8')
                # 创建sheet
                data_sheet = workbook.add_sheet('demo')

                # 列表格式数据
                # 定义循环下标
                index = 0
                for i in d:
                    # 每一列的内容(i)
                    for x, item in enumerate(i):
                        # 下标(x)，单元元素(item)
                        data_sheet.write(index, x, item, set_style('Times New Roman', 220, True))
                    index += 1
                # sys.exit();
                # 保存文件
                workbook.save('matrix_DJS.xls')

            if __name__ == '__main__':
                write_excel()
                print('创建matrix_DJS.xlsx文件成功')
# MakeMatrix_DJS().inputinfor()


class MakeMatrix_DJSa:
    def newyear(self):
        while True:
            h = int(input("h "))
            n = int(input("n "))
            a = float(input("a "))
            q = []
            d = [[] for i in range(n)]
            for i in range(n):
                q.append(prob(input("fname"), h))
            # print('@', q)

            for i in range(n):
                for j in range(n):
                    assert a > -1 and a < 1, 'result is nan'
                    result = 0
                    M = []
                    n0 = min(len(q[i]), len(q[j]))
                    for k in range(n0):
                        M.append((q[i][k] + q[j][k]) / 2)
                    # if n0 == len(L1):
                    #     for p in L2[n:]:
                    #         M.append(p / 2)
                    # elif n == len(L2):
                    #     for p in L1[n:]:
                    #         M.append(p / 2)
                    gam = math.gamma(a + 1)
                    # print(gam)
                    digam = digamma(1) - digamma(1 - a)
                    # print(digam)
                    for num in q[i]:
                        if num != 0.0:
                            result -= num * (num ** (-a)) * (math.log(num) + digam)
                            # print(result,'@')
                    for num in q[j]:
                        if num != 0.0:
                            result -= num * (num ** (-a)) * (math.log(num) + digam)
                            # print(result,'#')
                    result = result / 2
                    # print(result)
                    for num in M:
                        if num != 0.0:
                            result += num * ((num ** (-a))) * (math.log(num) + digam)
                    d[i].append(result / gam)

            # return d
            # print('@', d)
            #         print('@', result/2)
            #         print("--------------------")
            #     print("--------------------")
            # print("--------------------")
            def set_style(name, height, bold=False):
                style = xlwt.XFStyle()  # 初始化样式
                font = xlwt.Font()  # 为样式创建字体
                font.name = name
                font.bold = bold
                font.color_index = 4
                font.height = height
                style.font = font
                return style

            def write_excel():
                # 创建工作簿
                workbook = xlwt.Workbook(encoding='utf-8')
                # 创建sheet
                data_sheet = workbook.add_sheet('demo')

                # 列表格式数据
                # 定义循环下标
                index = 0
                for i in d:
                    # 每一列的内容(i)
                    for x, item in enumerate(i):
                        # 下标(x)，单元元素(item)
                        data_sheet.write(index, x, item, set_style('Times New Roman', 220, True))
                    index += 1
                # sys.exit();
                # 保存文件
                workbook.save('matrix_DJSa.xls')

            if __name__ == '__main__':
                write_excel()
                print('创建matrix_DJSa.xlsx文件成功')
# MakeMatrix_DJSa().newyear()

class MakeMatrix_DJSa:
    def newyear(self):
        while True:
            h = int(input("h "))
            n = int(input("n "))
            a = float(input("a "))
            q = []
            d = [[] for i in range(n)]
            for i in range(n):
                q.append(prob(input("fname"), h))
            # print('@', q)

            for i in range(n):
                for j in range(n):
                    assert a > -1 and a < 1, 'result is nan'
                    result = 0
                    M = []
                    n0 = min(len(q[i]), len(q[j]))
                    for k in range(n0):
                        M.append((q[i][k] + q[j][k]) / 2)
                    # if n0 == len(L1):
                    #     for p in L2[n:]:
                    #         M.append(p / 2)
                    # elif n == len(L2):
                    #     for p in L1[n:]:
                    #         M.append(p / 2)
                    gam = math.gamma(a + 1)
                    # print(gam)
                    digam = digamma(1) - digamma(1 - a)
                    # print(digam)
                    for num in q[i]:
                        if num != 0.0:
                            result -= num * (num ** (-a)) * (math.log(num) + digam)
                            # print(result,'@')
                    for num in q[j]:
                        if num != 0.0:
                            result -= num * (num ** (-a)) * (math.log(num) + digam)
                            # print(result,'#')
                    result = result / 2
                    # print(result)
                    for num in M:
                        if num != 0.0:
                            result += num * ((num ** (-a))) * (math.log(num) + digam)
                    d[i].append(-result / gam)

            # return d
            # print('@', d)
            #         print('@', result/2)
            #         print("--------------------")
            #     print("--------------------")
            # print("--------------------")
            def set_style(name, height, bold=False):
                style = xlwt.XFStyle()  # 初始化样式
                font = xlwt.Font()  # 为样式创建字体
                font.name = name
                font.bold = bold
                font.color_index = 4
                font.height = height
                style.font = font
                return style

            def write_excel():
                # 创建工作簿
                workbook = xlwt.Workbook(encoding='utf-8')
                # 创建sheet
                data_sheet = workbook.add_sheet('demo')

                # 列表格式数据
                # 定义循环下标
                index = 0
                for i in d:
                    # 每一列的内容(i)
                    for x, item in enumerate(i):
                        # 下标(x)，单元元素(item)
                        data_sheet.write(index, x, item, set_style('Times New Roman', 220, True))
                    index += 1
                # sys.exit();
                # 保存文件
                workbook.save('matrix_DJSa.xls')

            if __name__ == '__main__':
                write_excel()
                print('创建matrix_DJSa.xlsx文件成功')
# MakeMatrix_DJSa().newyear()


class MakeMatrix_DJS_PE:
    def newyear(self):
        while True:
            n = int(input("n "))
            N = int(input("N "))
            data = int(input("data "))
            q = []
            d = [[] for i in range(n)]
            for i in range(n):
                q.append(prob_pentropy(input("fname "), N, data))
            # print('@', q)

            for i in range(n):
                for j in range(n):
                    M = q[i]
                    for num in q[j].keys():
                        M[num] = (M[num] + q[j][num]) / 2
                    # print(M)
                    prob_d1 = list(filter(lambda x: x != 0, list(map(lambda k: q[i][k], q[i].keys()))))
                    prob_d2 = list(filter(lambda x: x != 0, list(map(lambda k: q[j][k], q[j].keys()))))
                    prob_m = list(filter(lambda x: x != 0, list(map(lambda k: M[k], M.keys()))))
                    pe_1 = -sum(prob_d1 * (numpy.log(prob_m) - numpy.log(prob_d1)))
                    pe_2 = -sum(prob_d2 * (numpy.log(prob_m) - numpy.log(prob_d2)))
                    # pe_m = -sum(prob_m * numpy.log(prob_m))
                    result = 1 / 2 * (pe_1 + pe_2)
                    d[i].append(result)

            # return d
            # print('@', d)
            #         print('@', result/2)
            #         print("--------------------")
            #     print("--------------------")
            # print("--------------------")
            def set_style(name, height, bold=False):
                style = xlwt.XFStyle()  # 初始化样式
                font = xlwt.Font()  # 为样式创建字体
                font.name = name
                font.bold = bold
                font.color_index = 4
                font.height = height
                style.font = font
                return style

            def write_excel():
                # 创建工作簿
                workbook = xlwt.Workbook(encoding='utf-8')
                # 创建sheet
                data_sheet = workbook.add_sheet('demo')

                # 列表格式数据
                # 定义循环下标
                index = 0
                for i in d:
                    # 每一列的内容(i)
                    for x, item in enumerate(i):
                        # 下标(x)，单元元素(item)
                        data_sheet.write(index, x, item, set_style('Times New Roman', 220, True))
                    index += 1
                # sys.exit();
                # 保存文件
                workbook.save('matrix_DJS_PE.xls')

            if __name__ == '__main__':
                write_excel()
                print('创建matrix_DJS_PE.xlsx文件成功')
# MakeMatrix_DJS_PE().newyear()

class MakeMatrix_DJSa_PE:
    def newyear2(self):
        while True:
            n = int(input("n "))
            N = int(input("N "))
            data = int(input("data "))
            a = float(input("a "))
            q = []
            d = [[] for i in range(n)]
            mi=[]
            ma=[]
            for i in range(n):
                q.append(prob_pentropy(input("fname "), N, data))
            # print('@', q)

            for i in range(n):
                for j in range(n):
                    assert a > -1 and a < 1, 'result is nan'
                    result = 0
                    M = []
                    n0 = min(len(q[i]), len(q[j]))
                    M = q[i]
                    for num in q[j].keys():
                        M[num] = (M[num] + q[j][num]) / 2
                    # if n0 == len(L1):
                    #     for p in L2[n:]:
                    #         M.append(p / 2)
                    # elif n == len(L2):
                    #     for p in L1[n:]:
                    #         M.append(p / 2)
                    prob_d1 = list(filter(lambda x: x != 0, list(map(lambda k: q[i][k], q[i].keys()))))
                    prob_d2 = list(filter(lambda x: x != 0, list(map(lambda k: q[j][k], q[j].keys()))))
                    prob_m = list(filter(lambda x: x != 0, list(map(lambda k: M[k], M.keys()))))
                    gam = math.gamma(a + 1)
                    # print(gam)
                    digam = digamma(1) - digamma(1 - a)
                    # print(digam)
                    for num in prob_d1:
                        if num != 0.0:
                            result -= num * (num ** (-a)) * (math.log(num) + digam)
                            # print(result,'@')
                    for num in prob_d2:
                        if num != 0.0:
                            result -= num * (num ** (-a)) * (math.log(num) + digam)
                            # print(result,'#')
                    result = result / 2
                    # print(result)
                    for num in prob_m:
                        if num != 0.0:
                            result += num * ((num ** (-a))) * (math.log(num) + digam)
                    d[i].append(abs(result / gam))
            for s in d:
                mi.append(min(s))
                ma.append(max(s))
            mifinal=min(mi)
            mafinal=max(ma)

            for s in d:
                for l in s:
                    l=(l-mifinal)/(mafinal-mifinal)

            # return d
            # print('@', d)
            #         print('@', result/2)
            #         print("--------------------")
            #     print("--------------------")
            # print("--------------------")
            def set_style(name, height, bold=False):
                style = xlwt.XFStyle()  # 初始化样式
                font = xlwt.Font()  # 为样式创建字体
                font.name = name
                font.bold = bold
                font.color_index = 4
                font.height = height
                style.font = font
                return style

            def write_excel():
                # 创建工作簿
                workbook = xlwt.Workbook(encoding='utf-8')
                # 创建sheet
                data_sheet = workbook.add_sheet('demo')

                # 列表格式数据
                # 定义循环下标
                index = 0
                for i in d:
                    # 每一列的内容(i)
                    for x, item in enumerate(i):
                        # 下标(x)，单元元素(item)
                        data_sheet.write(index, x, item, set_style('Times New Roman', 220, True))
                    index += 1
                # sys.exit();
                # 保存文件
                workbook.save('matrix_DJSa_PE.xls')

            if __name__ == '__main__':
                write_excel()
                print('创建matrix_DJSa_PE.xlsx文件成功')
# MakeMatrix_DJSa_PE().newyear2()
