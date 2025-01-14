# 2024/05

import matplotlib.pyplot as plt
import matplotlib.dates as mdates  
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import numpy as np

def draw_re(Ns,Es,Us,method=False):
    # 设置绘图风格  
    plt.style.use('seaborn-darkgrid')  
    
    # 创建一个新的图形  
    plt.figure(figsize=(12,6))  
    
    # 绘制x值  
    plt.plot(np.arange(len(Ns)), Ns, label='Ns', color='#A5C496', linewidth=2.0) 
    plt.plot(np.arange(len(Ns)), Es, label='Es', color='#C7988C', linewidth=2.0)
    plt.plot(np.arange(len(Ns)), Us, label='Us', color='#8891DB', linewidth=2.0)     

    # 设置标题和轴标签  
    # plt.title(f'精密星历与广播星历误差对比',fontproperties = 'SimSun',fontsize=20)  
    plt.xlabel('历元',fontproperties = 'SimSun')  
    plt.ylabel('坐标（米）',fontproperties = 'SimSun')  
    # 显示图例  
    plt.legend()  
    plt.xticks(rotation=45,fontproperties = 'Times New Roman',fontsize=15)  # 如果时间标签重叠，可以旋转它们以便更好地显示 
    plt.yticks(fontproperties = 'Times New Roman',fontsize=20)  
    # plt.autofmt_xdate()  # 自动调整日期标签的位置，避免重叠  
    
    # 显示图形  
    plt.tight_layout()
    if method:   
        plt.savefig(f'NEU_iter.pdf') 
    plt.savefig(f'NEU.pdf')
    plt.show()

def draw_pdf(Ms):
    # 设置绘图风格  
    plt.style.use('seaborn-darkgrid')  
    
    # 创建一个新的图形  
    plt.figure(figsize=(12,6))  
    
    # 绘制x值  
    plt.plot(np.arange(len(Ms)), Ms, color='#A5C496', linewidth=2.0)   

    # 设置标题和轴标签  
    # plt.title(f'精密星历与广播星历误差对比',fontproperties = 'SimSun',fontsize=20)  
    plt.xlabel('历元',fontproperties = 'SimSun')  
    plt.ylabel('坐标（米）',fontproperties = 'SimSun')  
    # 显示图例  
    plt.legend()  
    
    # 格式化x轴，使其显示合适的时间戳格式  
    plt.xticks(rotation=45,fontproperties = 'Times New Roman',fontsize=15)
    plt.yticks(fontproperties = 'Times New Roman',fontsize=20)  
    # plt.autofmt_xdate()  # 自动调整日期标签的位置，避免重叠  
    
    # 显示图形  
    plt.tight_layout()   
    plt.savefig(f'NEU.pdf') 
    plt.show()

def draw_NE(Ns,Es,Us,method=False):
    # 设置绘图风格  
    # plt.style.use('seaborn-darkgrid')  
    
    # 创建一个新的图形  
    plt.figure(figsize=(12,6))  
    plt.axes().set_aspect('equal', adjustable='box') 

    
    # 绘制x值  
    plt.scatter(Es,Ns, label='point', color='#A5C496', alpha=0.5,marker="+",s=20)
    plt.scatter(0,0,alpha=0.5,label='origin',s=25)

    # 设置标题和轴标签  
    # plt.title(f'精密星历与广播星历误差对比',fontproperties = 'SimSun',fontsize=20)  
    plt.xlabel('E',fontproperties = 'Times New Roman',fontsize=15)  
    plt.ylabel('N',fontproperties = 'Times New Roman',fontsize=15)  
    # 显示图例  
    plt.legend()  
    # 格式化x轴，使其显示合适的时间戳格式  
    # plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))    
    plt.xticks(fontproperties = 'Times New Roman',fontsize=15)  # 如果时间标签重叠，可以旋转它们以便更好地显示 
    plt.yticks(fontproperties = 'Times New Roman',fontsize=15)  
    # plt.autofmt_xdate()  # 自动调整日期标签的位置，避免重叠  
    
    # 显示图形  
    plt.tight_layout()   
    if method:
        plt.savefig(f'NE_scatter_iter.pdf') 

    plt.savefig(f'NE_scatter.pdf') 
    plt.show()