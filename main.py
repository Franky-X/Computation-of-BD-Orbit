import ephem
import re
import math
from datetime import datetime, timedelta
import numpy as np

# Define the GPS epoch (January 6, 1980)
gps_epoch = datetime(1980, 1, 6)

#RINEX P91 参数对应表
#北斗信号接口控制文件 P35
#https://zhuanlan.zhihu.com/p/111325516?from_voters_page=true

# 打开RINEX文件
rinex_file = "BRD400DLR_S_20230010000_01D_MN.rnx"
obs_data = open(rinex_file, "r")

# 设置计算时间
compute_time = ephem.Date("2023-01-01 01:00:00")

# 定义一个函数，用于提取特定时间的观测数据
def extract_observation_data(prn, time):
    observation_data = []

    for line in obs_data:
        if line.startswith(f"> EPH {prn} CNV1"):
            line = obs_data.readline()
            time_str = line.split()
            time_str = time_str[1] + "-" + time_str[2] + "-" + time_str[3] + " " + time_str[4] + ":"  + time_str[5] + ":" + time_str[6] 
            # print(time_str)
            obs_time = ephem.Date(time_str)
            if obs_time == time:
                for i in range(10):
                    observation_data.append(line)
                    line = obs_data.readline()
    
    return observation_data

# 指定要查找的卫星编号（北斗导航卫星C30的卫星编号）
satellite_prn = "C30"

# 提取特定时间的观测数据
observation_data = extract_observation_data(satellite_prn, compute_time)
broadcast_orbit = []

if observation_data:
    # 打印观测数据
    print(f"北斗导航卫星{satellite_prn}在{compute_time}的观测数据:")
    for line in observation_data:
        print(line)
        pattern = r"-?\d+\.\d+e[+-]\d+"
        matches = re.findall(pattern, line)
        # 将匹配到的字符串转换为浮点数
        float_numbers = [float(match) for match in matches]
        broadcast_orbit.append(float_numbers)
        # print(float_numbers)
else:
    print(f"找不到卫星{satellite_prn}在{compute_time}的观测数据")

#计算轨道参数
#BDCS 坐标系下的地心引力常数 
u = 3.986004418e+14
#BDCS 坐标系下的地球自转角速度
omega_e_dot = 7.2921150e-05
#圆周率
pi = 3.1415926535898
#计算长半轴
A = broadcast_orbit[2][3] * broadcast_orbit[2][3]
print("A:",A)

#计算卫星平均角速度
n0 = math.sqrt(u/(A*A*A))
print("n0:",n0)
#计算观测历元到参考历元的时间差
desired_datetime = datetime(2023, 1, 1, 1, 0, 30) 
# Calculate the time difference from the GPS epoch
time_difference = desired_datetime - gps_epoch

# Calculate the number of weeks
gps_week = time_difference.days // 7

# Calculate the time of day (in seconds)
time_of_day_seconds = time_difference.seconds + time_difference.days * 24 * 60 * 60
week = time_of_day_seconds / 60 / 60 / 24 / 7
time_of_week = (week - gps_week) * 604800
toe = broadcast_orbit[3][0]
toc = broadcast_orbit[3][0]
tk = time_of_week - broadcast_orbit[3][0]
print("toc",time_of_week)
t = time_of_week

if (tk > 302400):
    tk = tk-604800
if (tk < -302400):
    tk = tk+604800
print("tk",tk)

#改正平均角速度
n = n0 + broadcast_orbit[1][2]
print("n:",n)

#计算平近点角
Mk = broadcast_orbit[1][3] - n*tk

#迭代计算偏近点角 此处迭代三次
e = broadcast_orbit[2][1]
Ek = Mk
for i in range(30):
    Ek = Ek + e * math.sin(Ek)
print("Mk:",Mk)
print("Ek:",Ek)

#计算真近点角
vk = math.atan2((math.sqrt(1-e*e))*math.sin(Ek)/(1-e*math.cos(Ek)), (math.cos(Ek)-e)/(1-e*math.cos(Ek)))
print("vk:",vk)

#计算纬度幅角
phi_k = vk + broadcast_orbit[4][2]
print("phi_k:",phi_k)

#纬度幅角改正项
#径向改正项
#轨道倾角改正项
Cus = broadcast_orbit[2][2]
Cuc = broadcast_orbit[2][0]
Crs = broadcast_orbit[1][1]
Crc = broadcast_orbit[4][1]
Cis = broadcast_orbit[3][3]
Cic = broadcast_orbit[3][1]

du = Cus*math.sin(2*phi_k) + Cuc*math.cos(2*phi_k)
dr = Crs*math.sin(2*phi_k) + Crc*math.cos(2*phi_k)
di = Cis*math.sin(2*phi_k) + Cic*math.cos(2*phi_k)
print("du:",du)
print("dr:",dr)
print("di:",di)

#计算改正后的纬度幅角
uk = phi_k + du
print("uk:",uk)


#计算改正后的径向
rk = A*(1-e*math.cos(Ek)) + dr
print("rk:",rk)


#计算改正后的轨道倾角
i0 = broadcast_orbit[4][0]
IDOT = broadcast_orbit[5][0]
ik =  i0 + IDOT*tk + di
print("ik:",ik)

#计算卫星在轨道平面内的坐标
xk = rk*math.cos(uk)
yk = rk*math.sin(uk)
print("xk:",xk)
print("yk:",yk)

#计算历元升交点经度（地固系）
#计算MEO/IGSO 卫星在BDCS 坐标系中的坐标
omega_0 = broadcast_orbit[3][2]
omega_dot = broadcast_orbit[4][3]
omega_k = omega_0 + (omega_dot-omega_e_dot)*tk - omega_e_dot*toe
print("omega_k",omega_k)

Xk = xk * math.cos(omega_k) - yk * math.cos(ik)*math.sin(omega_k)
Yk = xk * math.sin(omega_k) + yk * math.cos(ik)*math.cos(omega_k)
Zk = yk * math.sin(ik)
print("BDCS:",Xk,Yk,Zk)

#计算历元升交点经度（惯性系）
#计算GEO 卫星在自定义坐标系中的坐标
#计算GEO 卫星转换至BDCS 坐标系中的坐标
omega_k = omega_0 + (omega_dot)*tk - omega_e_dot*toe
XGk = xk * math.cos(omega_k) - yk * math.cos(ik)*math.sin(omega_k)
YGk = xk * math.sin(omega_k) + yk * math.cos(ik)*math.cos(omega_k)
ZGk = yk * math.sin(ik)

phi = -5 / 360 * 2 * pi
RX = np.array([[1, 0, 0],
                  [0, math.cos(phi), math.sin(phi)],
                  [0, -math.sin(phi), math.cos(phi)]])
vector = np.array([XGk, YGk, ZGk])
result = np.dot(RX, vector)
phi = omega_e_dot*tk
RZ = np.array([[math.cos(phi), math.sin(phi), 0],
                  [-math.sin(phi), math.cos(phi), 0],
                  [0, 0, 1]])

result = np.dot(RZ, result)
print(result)

a0 = broadcast_orbit[0][0]
a1 = broadcast_orbit[0][1]
a2 = broadcast_orbit[0][2]

#计算卫星钟差
c = 3e8
F = -2 * math.sqrt(u) / (c*c) 
delta_tr = F * e * broadcast_orbit[2][3] * math.sin(Ek)
delta_tsv = a0 + a1*(t - toc) + a2*(t - toc)*(t - toc) + delta_tr
print("delta_tsv",delta_tsv)
# 关闭RINEX文件
obs_data.close()
