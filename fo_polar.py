# 注意要将mean值缩短到小数点后5位，否则会溢出

#coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def result_pic(result):
    """
    雷达图的绘制
    :param result: 分类数据
    :return: 雷达图
    """
    # 解析出类别标签和种类
    labels = ['PHG','VC','PCC','mCC','Hipp','Insula','MFG','MPFC','IFG','ACC','Crus','Ver','Md','Th','St','SM','PMC','STG','MTG','ITC']
    #labels = ['PCC', 'vmPFC', 'AG_R', 'AG_L', 'dlPFC_R', 'dlPFC_L', 'IPS_R', 'IPS_L',
           #   'dACC', 'insula_R', 'insula_L']
    kinds = list(result.iloc[:, 0])

    # 由于在雷达图中，要保证数据闭合，这里就再添加L列，并转换为 np.ndarray
    result = pd.concat([result, result[['PHG']]], axis=1)
   # result = pd.concat([result, result[['insula_L']]], axis=1)
    centers = np.array(result.iloc[:, 1:])

    # 分割圆周长，并让其闭合
    n = len(labels)
    angle = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angle = np.concatenate((angle, [angle[0]]))

    # 绘图
    
    # plt.rcParams["grid.linestyle"] = (5,9)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar') # 参数polar, 以极坐标的形式绘制图形
    
    
    # 画线
    ax.plot(angle, centers[0], linewidth=1,color=(65/256,109/256,180/256),label=kinds[0])
    ax.plot(angle, centers[1], linewidth=1,color=(180/256,29/256,35/256),label=kinds[1])
      #  lc = LineCollection(angle, linewidths=2, color=color)
      #  ax.add_collection(lc)
        # ax.fill(angle, centers[i])  # 填充底色

    # 添加属性标签
    plt.grid(b=True, which='both', color='0.65',linewidth = '0.3',linestyle='--')
    ax.set_thetagrids(angle * 180 / np.pi,labels)
    plt.title('')
    plt.legend(loc='lower right')
    plt.ylim(-0.1,0.9)
    # plt.savefig('E://state.eps')
    plt.show()

    

if __name__ == '__main__':
    result = pd.read_csv('/picture_MGZ_segregation.csv', sep=',')
    result_pic(result)


