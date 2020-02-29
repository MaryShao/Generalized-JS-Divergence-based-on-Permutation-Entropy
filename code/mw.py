import math
import numpy
import scipy
import itertools
from scipy.special import digamma
import matplotlib.pyplot as plt
#import xlwt


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
    # global x1, y1
    # x1=0
    # y1=0
    for l in final_list:
        # print(l)
        x = float(l[5])
        y = float(l[6])
        if not(math.isnan(x)) and not(math.isnan(y)) and y!=0.0 and x!=0.0:
        #     x1 = x
        #     y1 = y
        # print(x1,y1)
        #     print(x,y)
            point_l.append(Point(float(x), float(y)))
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

def prob(fname,h):   #square
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
    field = max(max_x,max_y)
    r=field/h
    while k-1<= h:
        n = 0
        for point in point_l:
            if (-k * r) <= point._x <= k * r and (-k * r) <= point._y <= k * r:
                n += 1
        p = (n - n_last) / N
        L1.append(p)
        k += 1
        n_last = n
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
    if l != []:
        list_x = []
        list_y = []
        for point in l:
            list_x.append(point._x)
            list_y.append(point._y)
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

    
def DJS(d1,d2,h):
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
    return x  #y

    #print(x)
    #print(y)
'''
x = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100,4200]
y1 = [0.2495862344768682, 0.6862157087543463, 0.22346540475411977, 0.37334048520872465, 0.4942936418571346, 0.3346175391914093, 0.20191797817578067, 0.18588695317592507, 0.20820379241889697, 0.24108214560217667, 0.19797337505395732, 0.24623461156305854, 0.4680506470713249, 0.3213774746821013, 0.30525333542637306, 0.22726025819897466, 0.1966907983522364, 0.4453226957762215, 0.4929407588343903, 0.24317322832245825, 0.30021243130034975, 0.2136008858608507, 0.36803875599981645, 0.30236602317664685, 0.662458945273515, 0.22567282966587915, 0.6446268779207496, 0.2824041131406706, 0.2409221146486495, 0.3574186909701193, 0.3090194600156836, 0.18042241033947917, 0.4626881562198354, 0.32612528220838427, 0.5185415328668606, 0.22194546503949708, 0.4040537275438408, 0.24332833047686458, 0.2891395444897842, 0.0, 0.0, 0.0,0.0]
y2 = [0.2141475082467662, 0.2539061950569278, 0.21739511632155062, 0.3538297050434078, 0.578419528112473, 0.6631511599639348, 0.13100631029743276, 0.4938019422074206, 0.21247303849855576, 0.23380513857218596, 0.33361594248217064, 0.24597554543296357, 0.606145415334871, 0.29250872932756705, 0.38592689033042077, 0.22636695634703102, 0.45873384394227185, 0.46243859993964276, 0.6341619150454587, 0.23430520383016493, 0.2251955686651883, 0.24297442596664673, 0.374614610610304, 0.19538868069322943, 0.679284236948747, 0.2746359785817845, 0.5756027723529447, 0.3104839048732862, 0.13678350213026125, 0.4401052251084175, 0.17960056681199035, 0.25798745275417173, 0.4086763060128986, 0.2546132212113507, 0.3680627044987415, 0.2198753602776474, 0.5471846704401274, 0.20529673990604627, 0.13017708868781164, 0.21445198726316048, 0.4305329919642515, 0.0,0.0]
y3=[0.4244638477733192, 0.3832205521641139, 0.2897934534321765, 0.36598414843675064, 0.6502149765535068, 0.6862157087543459, 0.22514958365807033, 0.5617953235987657, 0.2891417228289355, 0.5004188530259049, 0.4164695333123468, 0.5041074462500025, 0.6193722112334887, 0.4388973622693586, 0.5195237879246514, 0.3871896759056808, 0.48353806649365344, 0.49129726332345613, 0.6697365244243244, 0.2420142074170487, 0.3694376240296786, 0.2414490315999801, 0.3910787841704842, 0.2818837703272884, 0.6710894074470686, 0.5330311143407386, 0.6459797609434935, 0.5742183429721057, 0.43588597296583903, 0.32632885348145935, 0.5521773047533199, 0.21627792305079466, 0.5382736777491046, 0.35816284641689144, 0.5356813855292071, 0.21294935456344785, 0.17886297406381094, 0.2938036333817404, 0.2306775353593986, 0.34440725640474995, 0.3495144739566311, 0.14744852550914817, 0.22748031419696646]
y4=[]
y5=[]
for i in range(0,4300, 100):
    y4.append(0)
    y5.append(1)
    #x = np.array(x)
    #y = np.array(y)
plt.figure(figsize=(15,4))
plt.plot(x, y1,  marker='o',mfc='w', ms=6,label=u'HSI vs. SSE')
plt.plot(x, y2, marker='*',mfc='w', ms=6,label=u'HSI vs. N225')
plt.plot(x, y3, marker='D',mfc='w', ms=6,label=u'HSI vs. NDX')
plt.plot(x,y4,'--')
plt.plot(x,y5,'--')
plt.ylim(-0.1, 1.1)
plt.ylabel('window_DJS')
plt.xlabel('window')
plt.title("Similarity test using JSD") #图标题
plt.legend(loc='center left', bbox_to_anchor=(1,0.8))
plt.tight_layout()
plt.show()    
'''
 
    
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
        y.append(-result / gam) 
    return y
    # print(x)
    # print(y)
'''
# window_DJSa('^N225.csv','^HSI.csv',100,0.5,100)
x = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700]
y1=[4.619250095052837, 7.861450676486353, 4.46978683912596, 6.396651409160847, 8.036922385170781, 5.689317207863742, 4.5519043330406745, 3.824780810091838, 3.8778276516927286, 4.838204113448078, 4.063883073532974, 4.665330255936426, 7.161858184547411, 5.855661129316897, 6.247721469886164, 4.339576663125848, 3.1463723191276447, 8.318820505238502, 8.06952331462092, 4.907192533119179, 5.871270206134405, 3.7816513895109476, 6.473373782109313, 5.419927804018335, 10.411706604871956, 4.607266851713762, 10.81478132773926, 5.1738522142226895, 4.384795666206729, 6.388812054532058, 5.826687940293191, 3.261778431395999, 7.507163912327152, 5.901065143590166, 8.802034784175882, 4.871937648860557, 7.4186653597173295, 4.430555614560502]
y2=[4.292368794412478, 6.5242273356141, 9.116029950432713, 9.198596598018131, 2.8374425530256806, 8.376425689662469, 4.7341563129921935, 5.279292403109683, 5.428555996819297, 4.864342636988538, 10.127479735371695, 6.1186899410727165, 6.5459894710434146, 4.874751369894051, 7.635621293835352, 7.2988498463946785, 10.413351789519576, 4.399694364919647, 3.787280269654585, 4.143324609824845, 5.292873970894612, 3.6629906625818, 10.738645756095872, 5.484471356309004, 8.099105989117897, 5.197055279009888, 2.9145704592599966, 7.661191347746305, 2.6755081126368063, 4.437476359426021, 6.145518306165718, 4.576896339421285, 6.171385141059292, 4.121748900746891, 9.224497563470432, 3.89417505673402, 2.849189486012258, 3.981970880517729]
y3=[5.230388977728751, 9.19870456623854, 7.067009046430916, 8.644261350922783, 10.221345332990756, 8.311689480294943, 9.32940834990838, 7.546891377194894, 7.669848282608633, 7.909837223743332, 9.894105839828512, 4.636505039305997, 6.264248512551162, 4.444786678792921, 6.126812533262829, 4.636903634290302, 9.844823606074817, 9.72711625745579, 8.623540060261824, 9.547913771931496, 6.884171171622444, 5.455058750323358, 8.358071710186387, 4.1197938654880675, 7.775969849838826, 5.7595547796164235, 8.599311797470516, 4.200030615321347, 2.733502499978161, 4.495585788572728, 3.5214781515207934, 5.560226396622984, 5.710665408732009, 2.5304998532651712, 3.8635934633528257, 4.61708122060885, 5.329831843052362, 6.001694729274667]

plt.figure(figsize=(15,4))
plt.plot(x, y1,  marker='o',mfc='w', ms=6,label=u'HSI vs. SSE')
plt.plot(x, y2, marker='*',mfc='w', ms=6,label=u'HSI vs. N225')
plt.plot(x, y3, marker='D',mfc='w', ms=6,label=u'HSI vs. NDX')
plt.ylim(2, 12)
plt.ylabel('window_DJSa')
plt.xlabel('window')
plt.title("Similarity test using JSDa") #图标题
plt.legend(loc='center left', bbox_to_anchor=(1,0.8))
plt.tight_layout()
plt.show()      
'''  

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
    if point_l != []:
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
        if d1 != {} and d2 != {}:
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
    return Y

'''
x=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800]
y1=[0.003475135995240552, 0.004814448136127667, 0.00929170659320316, 0.011850410837854701, 0.0019709253342244437, 0.0024631521389654066, 0.0025077385250226146, 0.0019331692925305746, 0.001851130471569094, 0.0024652365332623744, 0.005890015478009625, 0.002302737330369552, 0.002469749561599288, 0.01062664248502896, 0.0029576642669872514, 0.008514644468547342, 0.001524716080888824, 0.008167402530116205, 0.0016217839119381959, 0.0002824562285061172, 0.01904725535370793, 0.008193757847866533, 0.005555022413628395, 0.0008901812744566018, 0.0013171791730135617, 0.003621190365786422, 0.003690162244428927, 0.0014799445566500355, 0.005177100680281233, 0.003766940191275974, 0.0024506006012366384, 0.003995847736844969, 0.007307960645525222, 0.008888273025909617, 0.0066296905242227856, 0.002742220650670965, 0.0006531070506358679, 0.0019534601691768433, 0.0008170017130155437]
y2=[0.017861002523314788, 0.003608678234492138, 0.006686587670000305, 0.0063809344329921296, 0.0005541223139330018, 0.005951908038820545, 0.013052995768413843, 0.0035136080454609636, 0.001371928971962838, 0.020963147885493742, 0.002457850073527175, 0.007273606325408647, 0.001276121313921259, 0.004525566643875965, 0.0026964588609810114, 0.0032696270409007916, 0.0064553594513676515, 0.006852819793480219, 0.0005146374688834105, 0.001541113865750272, 0.00437662761442752, 0.0007741910992179837, 0.0017540793293804487, 0.0024498855578874573, 0.0009464159276079496, 0.0012876149321239764, 0.0035180946653614606, 0.002677859643914353, 0.0033126246299399875, 0.00832434445843789, 0.0012814187954531147, 0.006767034531318418, 0.0027624127981824374, 0.0031050390646740443, 0.0006659327722868401, 0.005089190344125569, 0.0012970406396450905, 0.004566101809355981, 0.004471405499078106]
y3=[0.012651402183300055, 0.004967725414931219, 0.004830333877331883, 0.0057787525452805775, 0.004957240466096018, 0.0071457402986317, 0.001627733351605977, 0.0023536449914971754, 0.0014879083445871362, 0.0024450272949247116, 0.0017262574723817866, 0.005070844213466305, 0.0026003749679190374, 0.0023002632649045188, 0.004616992915249585, 0.010126308941318237, 0.0010571366441697954, 0.003324482090871632, 0.001490934174909555, 0.0024227384061368522, 0.002061677456982924, 0.003379560734183247, 0.006050208927696993, 0.003367817632511709, 0.005778689465109674, 0.0017854714819997416, 0.0014361811361410229, 0.008856405967406547, 0.0056367208148865346, 0.0023014002051457566, 0.00418215255311285, 0.0012131848035161815, 0.0009313895738296447, 0.00612994018203725, 0.006656068247560515, 0.012506088535644904, 0.015644937930915044, 0.0007691468954957822, 0.003493316221578103]

plt.figure(figsize=(15,4))
plt.plot(x, y1,  marker='o',mfc='w', ms=6,label=u'HSI vs. SSE')
plt.plot(x, y2, marker='*',mfc='w', ms=6,label=u'HSI vs. N225')
plt.plot(x, y3, marker='D',mfc='w', ms=6,label=u'HSI vs. NDX')
plt.ylim(0, 0.03)
plt.ylabel('window_PE_JSD')
plt.xlabel('window')
plt.title("Similarity test using PE_JSD") #图标题
plt.legend(loc='center left', bbox_to_anchor=(1,0.8))
plt.tight_layout()
plt.show() 
'''





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
        y.append(-result / gam) 
    return y

def window_DJSa_PE(fname_1,fname_2,a,n,data,lg):
    assert a > -1 and a < 1, 'result is nan'
    X = []
    Y = []
    gam = math.gamma(a+1)
    #print(gam)
    digam = digamma(1) - digamma(1 - a)    
    if isinstance(fname_1,list):
        l = max(len(fname_1),len(fname_2))
    else:
        l = max(len(in_point(fname_1)),len(in_point(fname_2)))
    for k in range(0,(l-lg), lg):
        M = []    
        d1 = window_prob_pentropy(fname_1, k, n, data,lg)
        d2 = window_prob_pentropy(fname_2, k, n, data,lg)
        M = d1
        if d1 != {} and d2 != {}:
            result=0
            for num in d2.keys():
                M[num] = (M[num] + d2[num]) / 2
            # print(M)
            # for i in range(0, len(M) - k): #not sure for moving window!
            prob_d1 = list(filter(lambda x: x!=0, list(map(lambda k: d1[k],d1.keys()))))
            prob_d2 = list(filter(lambda x: x!=0, list(map(lambda k: d2[k],d2.keys()))))
            prob_m = list(filter(lambda x: x!=0, list(map(lambda k: M[k], M.keys()))))
            pe_1 = -sum(prob_d1 * (numpy.log(prob_m) - numpy.log(prob_d1)))
            pe_2 = -sum(prob_d2 * (numpy.log(prob_m) - numpy.log(prob_d2)))
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
            X.append(k)
            Y.append(result/gam)
    return Y

x=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800,3900,4000]
y1=[-0.0034683073265366236, 0.002946699062046283, -0.0029254081412762624, -3.213651448022382e-05, 0.0056931774747574805, -0.00012728352566840462, 0.00015008505129218198, -0.0007427653497621262, -0.0008336744434153141, -0.0048774507252335635, 0.0003796144999613743, -0.0006643667812413127, -0.0009217418889283406, 0.0011304155020991735, -0.0014229613515485066, 0.00027943885246474275, -0.0006181916559794049, -6.915586438996928e-05, 0.00015359924204098933, 0.000559594051599661, -0.002617097613617851, -0.0009450512916730118, 0.0018284136220955772, -0.0003366166081768935, 0.0006289032779532384, 0.0017266397358166605, -0.003371798554672864, -0.001602375685613017, 0.000751560321242976, -0.0013611585941530693, -0.001899391496181906, -0.00023470161256119722, -2.436014588786713e-06, 0.005127633900462222, -0.002272118311824032, -5.5706215253193435e-05, -0.0002458321536003853, 0.0007753635743000014, 0.0014995423755664906,0.000001,-0.003]
y2=[-0.008529401205139078, 0.0008103658051358081, -0.0041899316786853285, 9.274180537181547e-05, -0.00011420393492888285, -0.0009280854748218488, -0.0034950997003624883, -0.0009925518208546542, -0.00026823886685287595, -0.00871405308692718, -0.0008788943099121332, -0.0003259180898539918, 0.0009802474129580454, -0.0009135960145504328, -0.002532261878468026, -0.001162464097006697, 0.00014328022478530758, -0.00017604447167761557, 6.069782003033638e-05, -8.958167447767637e-05, -0.0008925915446707137, -0.0012674756192358886, -0.0013613664781091522, -0.0005775984225902216, -7.91722101341398e-06, 0.00038146490524749587, -0.003225706432327113, -0.0011227238032280008, 0.0009790045809165689, -0.0045795023034492345, -0.00037929868241813695, -0.0015158784503175495, -0.00044573497275739475, 0.0013151590105983787, -0.0006114045923284515, 0.0004065010786961793, 0.00020934244417042977, 0.0017434337525871525, 0.00023252663452163545, -0.0013572819981899718, -0.0049194362481398335]
y3=[-0.0039558648405077996, 0.005018613441481313, -0.0032275893141721574, -0.0006569626780775375, 0.0030462227362228734, -2.434568529010977e-05, -0.0014973457842400856, -0.0037856198451109344, 0.0013722642633122297, -0.007034495238096572, 0.0007103601890603342, -0.001387363741266051, -0.0017495660771990084, 0.0023526054265257078, -0.0009910582346235634, 0.0006300917205004026, -9.530097328693524e-05, 0.0017497729186184351, 0.0004208321427753307, 0.00048254019871173196, -0.0022697769650515333, 0.00020618998180803726, 0.0019306887810133933, 6.75248903438929e-05, 0.00030658627423559584, -0.0010589952040479353, -0.0033119337416453445, -0.003662062937077071, -9.043286548465426e-06, -0.0017883116428255976, -0.0012336083202676255, -0.0004082499201386965, -0.00021669030080440955, 0.0018556702624758436, -0.00047907132952172856, -0.00011262412116887156, 2.7646079570484834e-05, 0.0022941834592238044, -0.0005284781542782721, 0.0019344117102030904, -0.005768033212110719]
plt.figure(figsize=(15,4))
plt.plot(x, y1,  marker='o',mfc='w', ms=6,label=u'HSI vs. SSE')
plt.plot(x, y2, marker='*',mfc='w', ms=6,label=u'HSI vs. N225')
plt.plot(x, y3, marker='D',mfc='w', ms=6,label=u'HSI vs. NDX')
plt.ylim(-0.01, 0.01)
plt.ylabel('window_Fractional_PE_JSD')
plt.xlabel('window')
plt.title("Similarity test using Fractional_PE_JSD") #图标题
plt.legend(loc='center left', bbox_to_anchor=(1,0.8))
plt.tight_layout()
plt.show()
