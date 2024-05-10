import numpy as np  
import math
from draw_picture import *

threshold=0.001
max_iterations=80

def read_file(filename):
    epochs = []  # 存储历元数据（每个历元是一个NumPy数组）
    times=[]  
    time=[]
    current_epoch_lines = []  # 存储当前历元的行数据
    with open(filename, 'r') as file: 
        for line in file:  
            line = line.strip()  # 去除行尾的换行符和空格  
            if line.startswith('#'):
                # time=line.split()[1:]  
                times.append(line)           
            if not line or line.startswith('#'):  # 如果行是空的或以#开头，则跳过  
                if current_epoch_lines:  # 如果当前历元有数据  
                    # 将当前历元的行数据转换为NumPy数组 
                    current_epoch = np.array([  
                        [float(field) for field in line.split()[1:]]  
                        for line in current_epoch_lines  
                    ])  
                    epochs.append(current_epoch)  # 将当前历元添加到epochs列表中  
                    current_epoch_lines = []  # 重置current_epoch_lines来收集下一个历元的行数据  
                continue  
            # 将当前行的数据添加到current_epoch_lines中  
            current_epoch_lines.append(line) 
        if current_epoch_lines:
             current_epoch = np.array([  
                        [float(field) for field in line.split()[1:]]  
                        for line in current_epoch_lines  
                    ])  
             epochs.append(current_epoch)  # 将当前历元添加到epochs列表中  
             current_epoch_lines = [] 

    return epochs,times
  
def cal_epoch(data,x0,y0,z0):
    t=0
    xi=data.T[0]
    yi=data.T[1]
    zi=data.T[2]
    ri=data.T[3]
    pi=data.T[4]
    P=np.diag(1/pi)
    delta=np.full((4,),np.inf) 
    iterations=0
    while np.linalg.norm(delta) > threshold and iterations < max_iterations:  
        r_calc = np.sqrt((xi-x0)**2 + (yi-y0)**2 + (zi-z0)**2)
        J = np.vstack([  
            (x0-xi) / r_calc,  
            (y0-yi) / r_calc,  
            (z0-zi) / r_calc,  
            np.ones_like(xi)
        ]).T  
        # 计算残差  
        residual = ri - (r_calc + t)  
        # 更新参数  
        # delta = np.linalg.lstsq(J, residual, rcond=None)[0]  
        delta = np.linalg.inv(J.T@ P @ J) @ J.T@ P @ residual
        x0 += delta[0]  
        y0 += delta[1]  
        z0 += delta[2]  
        t  += delta[3]
        iterations += 1 
    # test=(x**2+y**2)/6378140**2+z**2/6356755**2
    # print(x,y,z,test)
    rs=J @ delta.T-residual
    sigma= rs.T @ P @ rs
    sigma/=float(len(xi)-4)
    # print(len(xi))
    re=np.linalg.inv(J.T@ P @ J)*sigma
    return x0,y0,z0,t,sigma,re,len(xi)


def rad2degree(rad):
    return rad * 180 / np.pi

def xyz2BLH(x, y, z, rad=True):
    a=6378137.0000
    b=6356752.3142
    e2 = 1 - (b / a)**2
    p = np.sqrt(x**2+y**2)
    theta = np.arctan(z * a/(p * b))
    L = np.arctan2(y, x)
    B = np.arctan((z + e2/(1-e2)*b*np.sin(theta)**3)/(p - e2*a*np.cos(theta)**3))
    N = a/np.sqrt(1-e2*np.sin(B)**2)
    H = p / np.cos(B) - N
    if rad:
        return L, B, H
    else:
        return rad2degree(L), rad2degree(B),H


def cal_iter_all(epochs):
     print("start")
     xs=[]
     ys=[]
     zs=[]
     ts=[]
     x0,y0,z0,t,sigma,re,l=cal_epoch(epochs[0],0,0,0)# 定义初始x0 y0 z0
     xs.append(x0)
     ys.append(y0)
     zs.append(z0)
     ts.append(t)
     with open("result_iter_1.txt","w",encoding='utf-8') as f:
         f.write(f"{times[0]}\n")
         f.write(f"X(m)\t\t\tY(m)\t\t\tZ(m)\t\t\tT(m)\n")
         f.write(f"{x0:.4f}\t{y0:.4f}\t{z0:.4f}\t{t:.4f}\n")
         f.write(f"验后单位权中误差：{math.sqrt(sigma):.4f}(m)\n")
         f.write(f"验后估计方阵(m^2)\n")
            # f.write(f"{re}\n")
         for line in re:
              f_line = "\t".join(f"{item:.4e}" for item in line)  
              f.write(f"{f_line}\n")
         f.write("-------------------------------------------------------\n")

     index=1
    #  cal_iter(x0,y0,z0,t,sigma,re,l,epochs[1])
     for epoch in epochs[1:]:
          x1,y1,z1,t1,sigma1,re1,l1=cal_iter(x0,y0,z0,t,sigma,re,l,epoch)
          with open("result_iter_1.txt","a",encoding='utf-8') as f:
            f.write(f"{times[index]}\n")
            f.write(f"X(m)\t\t\tY(m)\t\t\tZ(m)\t\t\tT(m)\n")
            f.write(f"{x1:.4f}\t{y1:.4f}\t{z1:.4f}\t{t1:.4f}\n")
            f.write(f"验后单位权中误差：{math.sqrt(sigma1):.4f}(m)\n")
            f.write(f"验后估计方阵(m^2)\n")
                # f.write(f"{re}\n")
            for line in re1:
                f_line = "\t".join(f"{item:.4e}" for item in line)  
                f.write(f"{f_line}\n")
            f.write("-------------------------------------------------------\n")
          index+=1
          xs.append(x1)
          ys.append(y1)
          zs.append(z1)
          ts.append(t1)
          x0,y0,z0,t,sigma,re,l=x1,y1,z1,t1,sigma1,re1,l1

     xs=np.array(xs)
     ys=np.array(ys)
     zs=np.array(zs)
     ts=np.array(ts)
     print("end")
     return xs,ys,zs,ts
        


         

def cal_iter(x0,y0,z0,t,sigma,re,l,data):
     xi=data.T[0]
     yi=data.T[1]
     zi=data.T[2]
     ri=data.T[3]
     pi=data.T[4]
     P=np.diag(1/pi)
     delta=np.full((4,),np.inf) 
     Q=re/sigma
     iterations=0
     while np.linalg.norm(delta) > threshold and iterations < 1:  
        r_calc = np.sqrt((xi-x0)**2 + (yi-y0)**2 + (zi-z0)**2)
        h1 = np.vstack([  
                    (x0-xi) / r_calc,  
                    (y0-yi) / r_calc,  
                    (z0-zi) / r_calc,  
                    np.ones_like(xi)  
                ]).T 
        delta1=ri-(r_calc+t)
        K1=Q @ h1.T @ np.linalg.inv(np.linalg.inv(P)+h1@ Q @h1.T)
        delta= K1 @ delta1
        # 更新X  
        x0+=delta[0]
        y0+=delta[1]
        z0+=delta[2]
        t +=delta[3]

        Q1=Q-(K1 @ h1 @ Q)
        Q=Q1
        iterations+=1
     V1=h1@delta-delta1
     #  V1=h1@np.array([x0,y0,z0,t]).T-ri
     sigma1=(sigma*(l-4)+((delta1.T)@(K1.T)@np.linalg.inv(Q)@K1@delta1)+(V1.T@P@V1))/float(len(xi)-4)
    #  sigma1=(sigma*(l-4)+((delta1.T)@(K1.T)@Q@K1@delta1)+(V1.T@P@V1))/float(len(xi)-4)
     re1=Q1*sigma1
        
    #  print(iterations)
     return x0,y0,z0,t,sigma1,re1,len(xi)

def cal_all(epochs,times):
    print("start cal")
    x0=0
    y0=0
    z0=0
    index=0
    xs=[]
    ys=[]
    zs=[]
    ts=[]
    #清空文件
    with open('result.txt', 'w', encoding='utf-8') as f: 
        f.close()        
    for epoch in epochs:
        x,y,z,t,sigma,re,_=cal_epoch(epoch,x0,y0,z0)
        x0,y0,z0=x,y,z
        xs.append(x)
        ys.append(y)
        zs.append(z)
        ts.append(t)
        with open('result.txt', 'a', encoding='utf-8') as f: 
            f.write(f"{times[index]}\n")
            f.write(f"X(m)\t\t\tY(m)\t\t\tZ(m)\t\t\tT(m)\n")
            f.write(f"{x:.4f}\t{y:.4f}\t{z:.4f}\t{t:.4f}\n")
            f.write(f"验后单位权中误差：{math.sqrt(sigma):.4f}(m)\n")
            f.write(f"验后估计方阵(m^2)\n")
            # f.write(f"{re}\n")
            for line in re:
                f_line = "\t".join(f"{item:.4e}" for item in line)  
                f.write(f"{f_line}\n")
            f.write("-------------------------------------------------------\n")
        index+=1

    print("end cal")
    xs=np.array(xs)
    ys=np.array(ys)
    zs=np.array(zs)
    ts=np.array(ts)
    return xs,ys,zs,ts
        # print(x,y,z,t)

def cal_mean_std(xs,ys,zs,ts):
    Xmean=xs.mean()
    Ymean=ys.mean()
    Zmean=zs.mean()
    Tmean=ts.mean()
    Xstd=np.std(xs, ddof = 1)
    Ystd=np.std(ys, ddof = 1)
    Zstd=np.std(zs, ddof = 1)
    Tstd=np.std(ts, ddof = 1)
    # points=zip(xs,ys,zs)
    #输出mean std
    with open("mean.txt",'w') as meanf:
        meanf.write(f"{Xmean:.10e}&{Ymean:.10e}&{Zmean:.10e}&{Tmean:.10e}\n")
        meanf.write(f"{Xstd:.10e}&{Ystd:.10e}&{Zstd:.10e}&{Tstd:.10e}")

def cal_NEU(xs,ys,zs,x0,y0,z0):
    L,B,h=xyz2BLH(x0,y0,z0)
    Ns=[]
    Es=[]
    Us=[]
    #输出NEU
    with open("NEU.txt",'w') as xyz:
        for x,y,z in zip(xs,ys,zs):
            E=-(x-x0)*math.sin(L)+(y-y0)*math.cos(L)
            N=-(x-x0)*math.sin(B)*math.cos(L)-(y-y0)*math.sin(B)*math.sin(L)+(z-z0)*math.cos(B)
            U=(x-x0)*math.cos(B)*math.cos(L)+(y-y0)*math.cos(B)*math.sin(L)+(z-z0)*math.sin(B)
            xyz.write(f"{E:.4f}\t{N:.4f}\t{U:.4f}\n")
            Ns.append(N)
            Es.append(E)
            Us.append(U)
    Ns=np.array(Ns)
    Es=np.array(Es)
    Us=np.array(Us)
    return Ns,Es,Us

def cal_mean_RMS(Ns,Es,Us):
    Nmean=Ns.mean()
    Emean=Es.mean()
    Umean=Us.mean()
    Nrms=math.sqrt(sum(n**2 for n in Ns)/len(Ns))
    Erms=math.sqrt(sum(e**2 for e in Es)/len(Es))
    Urms=math.sqrt(sum(u**2 for u in Us)/len(Us))
    with open("NEU_iter.txt",'w') as f:
        f.write(f"{Nmean}\t{Emean}\t{Umean}\t")
        f.write(f"{Nrms}\t{Erms}\t{Urms}\t")
    # print(f"{Nmean}\t{Emean}\t{Umean}\t")
    # print(f"{Nrms}\t{Erms}\t{Urms}\t")


def cal(method=False):
    # xs,ys,zs,ts=cal_all(epochs,times)

    if method:
        xs,ys,zs,ts=cal_iter_all(epochs)
        cal_mean_std(xs,ys,zs,ts)
        Ns,Es,Us=cal_NEU(xs,ys,zs,x0,y0,z0)
        cal_mean_RMS(Ns,Es,Us)
        draw_re(Ns,Es,Us,method=True)
        draw_NE(Ns,Es,Us,method=True)
    else:
        xs,ys,zs,ts=cal_all(epochs,times)
        # cal_mean_std(xs,ys,zs,ts)
        Ns,Es,Us=cal_NEU(xs,ys,zs,x0,y0,z0)
        cal_mean_RMS(Ns,Es,Us)
        draw_re(Ns,Es,Us,method=True)
        draw_NE(Ns,Es,Us,method=True)


# 参考坐标 
x0=-1.13291501648681e6
y0=6.09252850388968e6
z0=1.50463316777129e6

epochs,times=read_file("CUSV_20212220_BDS_M0.5_I1.0_G2.0.txt")
cal()

# draw_pdf(Us)
# x0,y0,z0,t,sigma,re,l=cal_epoch(epochs[0],0,0,0)# 定义初始x0 y0 z0
# x1,y1,z1,t1,sigma1,re1,l1=cal_iter(x0,y0,z0,t,sigma,re,l,epochs[1])
# print(x1,y1,z1,t1,sigma1,l1)
# print(re1)
# cal_iter()
# xs,ys,zs,ts=cal_iter_all(epochs)
# print(ys.mean())
# print(np.std(ys, ddof = 1))


