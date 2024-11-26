# -*- coding: utf-8 -*-
"""
根据试验数据计算粘弹性参数
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import fmin

#广义Maxwell模型，离散本构
def analytical(paras, e, dt, dedt):
    Einf = paras[0]
    h=[]
    m = int(len(paras)/2)
    for i in range(m):
        h.append([0.])
        # print(h)
        Ei = paras[2*i+1]
        TAUi = paras[2*(i+1)]
        for j in range(1,len(dt)):
            h[i].append(np.exp(-dt[j]/TAUi)*h[i][j-1]+TAUi*Ei*(1-np.exp(-dt[j]/TAUi))*dedt[j])      
    ha=np.array(h).transpose()
    s = Einf*e+ha.sum(axis=1)
    return s
      
 #线性广义Maxwell模型，等应变速率解析解       
def analytical_eq(paras, e, t, de_eng_dt):
    Einf = paras[0]
    # alpha= paras[-2]
    # belta= paras[-1]
    h=0
    m = int(len(paras)/2)
    for i in range(m-1):
        Ei = paras[2*i+1]
        TAUi = paras[2*(i+1)]
        h+=de_eng_dt*Ei*TAUi*(1-np.exp(-t/TAUi))
    # s =  Einf*e+alpha*(e**2)+belta*(e**3)+h
    s =  Einf*e+h
    return s
          
#工程应力应变转化为真实应力应变
def cal_exp_data(e_eng, s_eng, de_eng_dt):
    e = np.log(e_eng + 1.0)
    s = s_eng * (1.0 + e_eng)
    t = e_eng / de_eng_dt
    e_temp = np.array([0])
    t_temp = np.array([0])
    e_temp = np.append(e_temp,e[:-1])
    t_temp = np.append(t_temp,t[:-1])
    dt = t - t_temp
    # de = e - e_temp
    dedt = 1/(t+1/de_eng_dt)
    e_sim =  np.log(de_eng_dt*t + 1.0)
    return e, s, t, dedt, dt, e_sim

#计算损失函数
def cal_cost(exp_data, paras):
    cost = 0
    for n in exp_data.keys():
        e_exp = exp_data[n]['e']
        s_exp = exp_data[n]['s']
        t = exp_data[n]['t'] 
        dedt = exp_data[n]['dedt']
        dt = exp_data[n]['dt']
        e_sim = exp_data[n]['e_sim']
        s_sim = analytical(paras, e_sim, dt, dedt)
        cost += np.sum((s_exp - s_sim)**2, axis=0)/len(s_exp)
    return cost

#目标函数=损失函数+罚函数
def func(x):
    cost = cal_cost(exp_data, x)
    punish = 0
    y = cost
    for i in range(len(x)-2):
        punish += (max(0, -x[i]))**2
    y += 10**100 * punish
    return y


def create_exp_data(gauge_len,disp_rate):
    exp_data = {}
    n = 0
    for disp_rate in disp_rate:
        data = pd.read_excel(io='time_sim.xlsx', sheet_name=str(disp_rate))#工程应力应变
        de_eng_dt = disp_rate/60.0/gauge_len

        for colomn in [['e1', 's1']]:
            n += 1

            e_eng = data[colomn[0]]
            s_eng = data[colomn[1]]

            e_eng = e_eng[e_eng.notna()]
            s_eng = s_eng[s_eng.notna()]
            
            max_index = np.argmax(s_eng)
            e_eng = e_eng[0:max_index]
            s_eng = s_eng[0:max_index]

            e_eng = np.array(e_eng - e_eng[0])
            s_eng = np.array(s_eng- s_eng[0])        

            e, s, t, dedt, dt, e_sim = cal_exp_data(e_eng, s_eng, de_eng_dt)

            exp_data[n] = {}
            exp_data[n]['e'] = e
            exp_data[n]['s'] = s
            exp_data[n]['t'] = t
            exp_data[n]['dedt'] = dedt
            exp_data[n]['dt'] = dt
            exp_data[n]['e_sim'] = e_sim
    return exp_data


def func(x):
    cost = cal_cost(exp_data, x)
    punish = 0
    y = cost
    for i in range(len(x)):
        punish += (max(0, -x[i]))**2
    y += 1e20 * punish
    return y


if __name__ == "__main__":
    gauge_len,disp_rate =70,[500]
    exp_data = create_exp_data(gauge_len,disp_rate)

    paras = [2,2,200,16,1]
    
    #####执行优化计算，
    # paras0 = [10,1,0.001,10,1,1,1000]#参数初始值
    # paras = fmin(func, paras0, maxiter=10000)

    #####参数的Prony级数形式
    E_inst = paras[0]+paras[1]+paras[3]
    gk1 = paras[1]/E_inst
    gk2 = paras[3]/E_inst
    # gk3 = paras[5]/E_inst
    prony=[E_inst,[gk1,paras[2]],[gk2,paras[4]]]

    print('umat',paras)
    print('prony',prony)
    #possion is 0.49
    for n in range(1,len(exp_data)+1):
        e_exp = exp_data[n]['e']
        s_exp = exp_data[n]['s']
        t = exp_data[n]['t']
        dedt = exp_data[n]['dedt']
        dt = exp_data[n]['dt']
        e_sim = exp_data[n]['e_sim']
        s_sim = analytical(paras, e_sim, dt, dedt)
        #输出数据
        np.savetxt("e.txt", e_exp)
        np.savetxt("s.txt", s_sim)
        np.savetxt("t.txt", t)
        # s_sim = analytical2(paras, dedt, t)
        if n <= 1:
            color='r'
        elif 2>=n>1:
            color='g'
        elif 3>=n>2:
            color='b'
        elif 4>=n>3:
            color='black'
        # plt.plot(e_exp, (s_exp-s_sim)/max(s_exp), color=color,linesteweyle=':')#归一化结果
        # plt.plot(e_exp, s_exp, color=color,linestyle=':')#试验值
        plt.plot(e_exp, s_sim,color=color,linestyle='solid',linewidth=2)#预测值
        plt.xlabel('true strain')
        plt.ylabel('true stress')

