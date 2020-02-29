import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from twod_DKL2 import DJS, DJSa, DJS_PE, DJSa_PE, in_point2, in_point
import numpy as np
from twod_DKL2 import DJS, DJSa, DJS_PE, DJSa_PE, in_point2, in_point
def stock_data(n, *vartuple):
    stock_in = []
    in_x=[]
    in_y=[]
    l=[]
    for var in vartuple:
        stock_in=in_point(var)
        l.append(len(stock_in))
        for point in stock_in:
            in_x.append(point._x)
            in_y.append(point._y)
    X = np.vstack((np.array(in_x),
                   np.array(in_y))).T
    # print(X)
    y = np.hstack([np.zeros(l[0], dtype=np.intp),
                   np.ones(l[1], dtype=np.intp),
                   np.ones(l[2], dtype=np.intp)*2,
                   np.ones(l[2], dtype=np.intp)*3,
                   np.ones(l[2], dtype=np.intp)*4,
                   np.ones(l[2], dtype=np.intp)*5,
                   np.ones(l[2], dtype=np.intp)*6,
                   np.ones(l[2], dtype=np.intp)*7,
                   np.ones(l[2], dtype=np.intp)*8])
                   #np.ones(l[2], dtype=np.intp)*9,
                   #np.ones(l[2], dtype=np.intp)*10,
                   #np.ones(l[2], dtype=np.intp)*11,
                   #np.ones(l[2], dtype=np.intp)*12,
                   #np.ones(l[2], dtype=np.intp)*13,
                   #np.ones(l[2], dtype=np.intp)*14])
                   #np.ones(l[2], dtype=np.intp)*15,
                   #np.ones(l[2], dtype=np.intp)*16,
                   #np.ones(l[2], dtype=np.intp)*17])
    return X, y

X, y = stock_data(9, "GS.csv","C.csv", "V.csv","JPM.csv","NMR.csv",
                  "PEP.csv","YUM.csv","KO.csv","JNJ.csv")
                  #"DIS.csv","NWS.csv",
                  #"GM.csv","F.csv","AEP.csv","CVX.csv")
                  #"NKE.csv","PG.csv","WMT.csv")

y=[ 0,  1,  2,  3,  4,  5,  6,  7,  8]


def jsd_kpca(X, N, k):
    # step 1
    K0=[]
    for i in range(N):
        for j in range(i):
            # print(Xl[i])
            K0.append(DJS(X[i], X[j]))
            # K0.append(DJSa(Xl[i], Xl[j], 100, 0.5))
            # K0.append(DJSa_PE(Xl[i], Xl[j], 0.5, 3))
            # K0.append(DJS_PE(Xl[i], Xl[j], 3))
            # K0.append((DJS_PE(Xl[i], Xl[j], 3)+DJS_PE(Xl[j], Xl[i], 3))/2)

    sq_dist=np.array(K0)
    K = squareform(sq_dist)
    # print(K)
    # step 2
    one_N = np.ones((N, N))/N
    K = K - one_N.dot(K) - K.dot(one_N) + one_N.dot(K).dot(one_N)
    # step 3
    Lambda, Q = np.linalg.eig(K)
    eigen_pairs = [(Lambda[i], Q[:, i]) for i in range(len(Lambda))]
    eigen_pairs = sorted(eigen_pairs, reverse=True, key=lambda k: k[0])
    # print(eigen_pairs)
    return np.column_stack((eigen_pairs[i][1] for i in range(k)))

def visual_kpca(X):
    #kpca_result
    X_kpca = jsd_kpca(X, 9, 2)  #k !effect??
    # print(X_kpca)

    #visual analysis
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    #
    ax[0].scatter(X_kpca[y==0, 0], X_kpca[y==0, 1],color='r', marker='o', alpha=.4,label='GS')
    ax[0].scatter(X_kpca[y==1, 0], X_kpca[y==1, 1],color='r', marker=',', alpha=.4,label='C')
    ax[0].scatter(X_kpca[y==2, 0], X_kpca[y==2, 1], color='r', marker='^', alpha=.4,label='V')
    ax[0].scatter(X_kpca[y==3, 0], X_kpca[y==3, 1], color='r', marker='+', alpha=.4,label='JPM')
    ax[0].scatter(X_kpca[y==4, 0], X_kpca[y==4, 1], color='r', marker='P', alpha=.4,label='NMR')
    ax[0].scatter(X_kpca[y==5, 0], X_kpca[y==5, 1], color='b', marker='*', alpha=.4,label='PEP')
    ax[0].scatter(X_kpca[y==6, 0], X_kpca[y==6, 1], color='b', marker='X', alpha=.4,label='YUM')
    ax[0].scatter(X_kpca[y==7, 0], X_kpca[y==7, 1], color='b', marker='D', alpha=.4,label='KO')
    ax[0].scatter(X_kpca[y==8, 0], X_kpca[y==8, 1], color='b', marker='4', alpha=.4,label='JNJ')
    #ax[0].scatter(X_kpca[y==9, 0], X_kpca[y==9, 1], color='y', marker='v', alpha=.4,label='DIS')
    #ax[0].scatter(X_kpca[y==10, 0], X_kpca[y==10, 1], color='y', marker='s', alpha=.4,label='NWS')
    #ax[0].scatter(X_kpca[y==11, 0], X_kpca[y==11, 1], color='g', marker='p', alpha=.4,label='GM')
    #ax[0].scatter(X_kpca[y==0, 0], X_kpca[y==12, 1],color='g', marker='1', alpha=.4,label='F')
    #ax[0].scatter(X_kpca[y==1, 0], X_kpca[y==13, 1],color='g', marker='2', alpha=.4,label='AEP')
    #ax[0].scatter(X_kpca[y==2, 0], X_kpca[y==14, 1], color='g', marker='3', alpha=.4,label='CVX')
    #ax[0].scatter(X_kpca[y==3, 0], X_kpca[y==15, 1], color='k', marker='4', alpha=.4,label='NKE')
    #ax[0].scatter(X_kpca[y==4, 0], X_kpca[y==16, 1], color='k', marker='d', alpha=.4,label='PG')
    #ax[0].scatter(X_kpca[y==5, 0], X_kpca[y==17, 1], color='k', marker='x', alpha=.4,label='WMT')
       

    label_count = np.bincount(y)
                                    # 统计各类别出现的次数
                                    # label_count[0] = 500
                                    # label_count[1] = 500
    ax[1].scatter(X_kpca[y==0, 0], np.zeros(label_count[0]),color='r', marker='o',label='GS')
    ax[1].scatter(X_kpca[y==1, 0], np.zeros(label_count[1]),color='r', marker=',', label='C')
    ax[1].scatter(X_kpca[y==2, 0], np.zeros(label_count[2]), color='r', marker='^',label='V')
    ax[1].scatter(X_kpca[y==3, 0], np.zeros(label_count[3]), color='r', marker='+', label='JPM')
    ax[1].scatter(X_kpca[y==4, 0], np.zeros(label_count[4]), color='r', marker='P', label='NMR')
    ax[1].scatter(X_kpca[y==5, 0], np.zeros(label_count[5]), color='b', marker='*', label='PEP')
    ax[1].scatter(X_kpca[y==6, 0], np.zeros(label_count[6]), color='b', marker='X', label='YUM')
    ax[1].scatter(X_kpca[y==7, 0], np.zeros(label_count[7]), color='b', marker='D', label='KO')
    ax[1].scatter(X_kpca[y==8, 0], np.zeros(label_count[8]), color='b', marker='4',label='JNJ')
    #ax[1].scatter(X_kpca[y==9, 0], np.zeros(label_count[0]), color='y', marker='v',label='DIS')
    #ax[1].scatter(X_kpca[y==10, 0], np.zeros(label_count[0]), color='y', marker='s',label='NWS')
    #ax[1].scatter(X_kpca[y==11, 0], np.zeros(label_count[0]), color='g', marker='p', label='GM')
    #ax[1].scatter(X_kpca[y==12, 0], np.zeros(label_count[0]),color='g', marker='1',label='F')
    #ax[1].scatter(X_kpca[y==13, 0], np.zeros(label_count[0]),color='g', marker='2', label='AEP')
    #ax[1].scatter(X_kpca[y==14, 0], np.zeros(label_count[0]), color='g', marker='3',label='CVX')
    #ax[1].scatter(X_kpca[y==15, 0], np.zeros(label_count[0]), color='k', marker='4', label='NKE')
    #ax[1].scatter(X_kpca[y==16, 0], np.zeros(label_count[0]), color='k', marker='d',label='PG')
    #ax[1].scatter(X_kpca[y==17, 0], np.zeros(label_count[0]), color='k', marker='x', label='WMT')
                                    # y轴置零
                                    # 投影到x轴
                                    
    ax[1].set_ylim([-1, 1])
    ax[0].set_xlabel('PC1')
    ax[0].set_ylabel('PC2')
    ax[1].set_xlabel('PC1')
    plt.legend(loc='center left', bbox_to_anchor=(1,0.8))
    plt.margins(0)
    plt.subplots_adjust(bottom=0.15)    
    plt.show()

# visual_pca(X)
print(jsd_kpca(X, 9, 2))
#visual_kpca(X)
# X_kpca = jsd_kpca(X, k=4)
# print(X_kpca[y==0, 0])
# print(X)