import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from oned_DKL2 import DJS, DJSa, DJS_PE, DJSa_PE, in_point2
def stock_data(*vartuple):
    stock_in = []
    in_x=[]
    for var in vartuple:
        stock_in.append(in_point2(var))
    l=min(len(data) for data in stock_in)
    for data in stock_in:
        in_x.append(data[0:l])

    X = (np.array([in_x[0],in_x[1],in_x[2],in_x[3],in_x[4],in_x[5],
                   in_x[6],in_x[7],in_x[8],in_x[9],in_x[10],in_x[11]]))
    # print(X)
    y = np.hstack([np.zeros(1, dtype=np.intp),
                   np.ones(1, dtype=np.intp),
                   np.ones(1, dtype=np.intp)*2,
                   np.ones(1, dtype=np.intp)*3,
                   np.ones(1, dtype=np.intp)*4,
                   np.ones(1, dtype=np.intp)*5,
                   np.ones(1, dtype=np.intp)*6,
                   np.ones(1, dtype=np.intp)*7,
                   np.ones(1, dtype=np.intp)*8,
                   np.ones(1, dtype=np.intp)*9,
                   np.ones(1, dtype=np.intp)*10,
                   np.ones(1, dtype=np.intp)*11])
    return X, y
X, y = stock_data("^HSI.csv","^DJI.csv", "^N225.csv",'^RUT.csv','^IXIC.csv',
                  '^GSPC.csv','^SSE.csv','^KS11.csv','^AEX.csv','^FCHI.csv',
                  '^GDAXI.csv','^BSESN.csv')
# print(X[y==2])

####################################################
##DJS_PE

def rbf_kpca(Xl, N, k):
    # step 1
    K0=[]
    for i in range(N):
        for j in range(i):
            print(Xl[i])
            # K0.append(DJS(Xl[i], Xl[j], 100))
            # K0.append(DJSa(Xl[i], Xl[j], 100, 0.5))
            K0.append(DJS_PE(Xl[i], Xl[j], 3))

    sq_dist=np.array(K0)
    K = squareform(sq_dist)
    print(K)
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
    X_kpca = rbf_kpca(X, len(X), 2)  #k !effect??
    # print(X_kpca)

    #visual analysis
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    #
    ax[0].scatter(X_kpca[y==0, 0], X_kpca[y==0, 1],color='r', marker='o', alpha=.4,label='^HSI')
    ax[0].scatter(X_kpca[y==1, 0], X_kpca[y==1, 1],color='b', marker=',', alpha=.4,label='^DJI')
    ax[0].scatter(X_kpca[y==2, 0], X_kpca[y==2, 1], color='r', marker='^', alpha=.4,label='^N225')
    ax[0].scatter(X_kpca[y==3, 0], X_kpca[y==3, 1], color='b', marker='+', alpha=.4,label='^RUT')
    ax[0].scatter(X_kpca[y==4, 0], X_kpca[y==4, 1], color='b', marker='P', alpha=.4,label='^IXIC')
    ax[0].scatter(X_kpca[y==5, 0], X_kpca[y==5, 1], color='b', marker='*', alpha=.4,label='^GSPC')
    ax[0].scatter(X_kpca[y==6, 0], X_kpca[y==6, 1], color='r', marker='X', alpha=.4,label='^SSE')
    ax[0].scatter(X_kpca[y==7, 0], X_kpca[y==7, 1], color='r', marker='D', alpha=.4,label='^KS11')
    ax[0].scatter(X_kpca[y==8, 0], X_kpca[y==8, 1], color='y', marker='4', alpha=.4,label='^AEX')
    ax[0].scatter(X_kpca[y==9, 0], X_kpca[y==9, 1], color='y', marker='v', alpha=.4,label='^FCHI')
    ax[0].scatter(X_kpca[y==10, 0], X_kpca[y==10, 1], color='y', marker='s', alpha=.4,label='^GDAXI')
    ax[0].scatter(X_kpca[y==11, 0], X_kpca[y==11, 1], color='r', marker='p', alpha=.4,label='^BSESN')

    label_count = np.bincount(y)
                                    # 统计各类别出现的次数
                                    # label_count[0] = 500
                                    # label_count[1] = 500
    ax[1].scatter(X_kpca[y==0, 0], np.zeros(label_count[0]),color='r', marker='o',label='^HSI')
    ax[1].scatter(X_kpca[y==1, 0], np.zeros(label_count[1]),color='b', marker=',',label='^DJI')
    ax[1].scatter(X_kpca[y==2, 0], np.zeros(label_count[2]), color='r', marker='^',label='^N225')
    ax[1].scatter(X_kpca[y==3, 0], np.zeros(label_count[3]), color='b', marker='+',label='^RUT')
    ax[1].scatter(X_kpca[y==4, 0], np.zeros(label_count[4]), color='b', marker='P',label='^IXIC')
    ax[1].scatter(X_kpca[y==5, 0], np.zeros(label_count[5]), color='b', marker='*',label='^GSPC')
    ax[1].scatter(X_kpca[y==6, 0], np.zeros(label_count[6]), color='r', marker='X',label='^SSE')
    ax[1].scatter(X_kpca[y==7, 0], np.zeros(label_count[7]), color='r', marker='D',label='^KS11')
    ax[1].scatter(X_kpca[y==8, 0], np.zeros(label_count[8]), color='y', marker='4',label='^AEX')
    ax[1].scatter(X_kpca[y==9, 0], np.zeros(label_count[9]), color='y', marker='v', label='^FCHI')
    ax[1].scatter(X_kpca[y==10, 0], np.zeros(label_count[10]), color='y', marker='s',label='^GDAXI')
    ax[1].scatter(X_kpca[y==11, 0], np.zeros(label_count[11]), color='r', marker='p',label='^BSESN')    
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
print(rbf_kpca(X, 3, 2))
# visual_kpca(Xl)
# X_kpca = rbf_kpca(X, k=4)
# print(X_kpca[y==0, 0])


####################################################


####################################################
##DJSa_PE

def rbf_kpca_aPE(Xl, N, k):
    # step 1
    K1=[]
    for i in range(N):
        for j in range(i):
            print(Xl[i])
            # K0.append(DJS(Xl[i], Xl[j], 100))
            # K0.append(DJSa(Xl[i], Xl[j], 100, 0.5))
            K1.append(DJSa_PE(Xl[i], Xl[j], 0.5,3))

    sq_dist=np.array(K1)
    K = squareform(sq_dist)
    print(K)
    # step 2
    one_N = np.ones((N, N))/N
    K = K - one_N.dot(K) - K.dot(one_N) + one_N.dot(K).dot(one_N)
    # step 3
    Lambda, Q = np.linalg.eig(K)
    eigen_pairs = [(Lambda[i], Q[:, i]) for i in range(len(Lambda))]
    eigen_pairs = sorted(eigen_pairs, reverse=True, key=lambda k: k[0])
    # print(eigen_pairs)
    return np.column_stack((eigen_pairs[i][1] for i in range(k)))

def visual_kpca_aPE(X):
    #kpca_result
    X_kpca = rbf_kpca_aPE(X, len(X), 2)  #k !effect??
    # print(X_kpca)

    #visual analysis
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    #
    ax[0].scatter(X_kpca[y==0, 0], X_kpca[y==0, 1],color='r', marker='o', alpha=.4,label='^HSI')
    ax[0].scatter(X_kpca[y==1, 0], X_kpca[y==1, 1],color='b', marker=',', alpha=.4,label='^DJI')
    ax[0].scatter(X_kpca[y==2, 0], X_kpca[y==2, 1], color='r', marker='^', alpha=.4,label='^N225')
    ax[0].scatter(X_kpca[y==3, 0], X_kpca[y==3, 1], color='b', marker='+', alpha=.4,label='^RUT')
    ax[0].scatter(X_kpca[y==4, 0], X_kpca[y==4, 1], color='b', marker='P', alpha=.4,label='^IXIC')
    ax[0].scatter(X_kpca[y==5, 0], X_kpca[y==5, 1], color='b', marker='*', alpha=.4,label='^GSPC')
    ax[0].scatter(X_kpca[y==6, 0], X_kpca[y==6, 1], color='r', marker='X', alpha=.4,label='^SSE')
    ax[0].scatter(X_kpca[y==7, 0], X_kpca[y==7, 1], color='r', marker='D', alpha=.4,label='^KS11')
    ax[0].scatter(X_kpca[y==8, 0], X_kpca[y==8, 1], color='y', marker='4', alpha=.4,label='^AEX')
    ax[0].scatter(X_kpca[y==9, 0], X_kpca[y==9, 1], color='y', marker='v', alpha=.4,label='^FCHI')
    ax[0].scatter(X_kpca[y==10, 0], X_kpca[y==10, 1], color='y', marker='s', alpha=.4,label='^GDAXI')
    ax[0].scatter(X_kpca[y==11, 0], X_kpca[y==11, 1], color='r', marker='p', alpha=.4,label='^BSESN')

    label_count = np.bincount(y)
                                    # 统计各类别出现的次数
                                    # label_count[0] = 500
                                    # label_count[1] = 500
    ax[1].scatter(X_kpca[y==0, 0], np.zeros(label_count[0]),color='r', marker='o',label='^HSI')
    ax[1].scatter(X_kpca[y==1, 0], np.zeros(label_count[1]),color='b', marker=',',label='^DJI')
    ax[1].scatter(X_kpca[y==2, 0], np.zeros(label_count[2]), color='r', marker='^',label='^N225')
    ax[1].scatter(X_kpca[y==3, 0], np.zeros(label_count[3]), color='b', marker='+',label='^RUT')
    ax[1].scatter(X_kpca[y==4, 0], np.zeros(label_count[4]), color='b', marker='P',label='^IXIC')
    ax[1].scatter(X_kpca[y==5, 0], np.zeros(label_count[5]), color='b', marker='*',label='^GSPC')
    ax[1].scatter(X_kpca[y==6, 0], np.zeros(label_count[6]), color='r', marker='X',label='^SSE')
    ax[1].scatter(X_kpca[y==7, 0], np.zeros(label_count[7]), color='r', marker='D',label='^KS11')
    ax[1].scatter(X_kpca[y==8, 0], np.zeros(label_count[8]), color='y', marker='4',label='^AEX')
    ax[1].scatter(X_kpca[y==9, 0], np.zeros(label_count[9]), color='y', marker='v', label='^FCHI')
    ax[1].scatter(X_kpca[y==10, 0], np.zeros(label_count[10]), color='y', marker='s',label='^GDAXI')
    ax[1].scatter(X_kpca[y==11, 0], np.zeros(label_count[11]), color='r', marker='p',label='^BSESN')    
                                    # y轴置零
                                    # 投影到x轴
                                    
    ax[1].set_ylim([-1, 1])
    ax[0].set_xlabel('PC1')
    ax[0].set_ylabel('PC2')
    ax[1].set_xlabel('PC1_aPE')
    plt.legend(loc='center left', bbox_to_anchor=(1,0.8))
    plt.margins(0)
    plt.subplots_adjust(bottom=0.15)    
    plt.show()
####################################################