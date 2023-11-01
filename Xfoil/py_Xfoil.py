import os.path
import subprocess
import numpy as np
from scipy.optimize import fsolve
import re

def read_pack(filename):
    '''
        Чтение данных из файла, записываемого Xfoil

        Ввод: filename (str) - имя файла, записанного Xfoil
        Вывод: float_data (numpy.ndarray) - массив данных
    '''

    file = open(filename)
    values = file.read().split("\n")
    data = []
    
    for key in values:
        value = re.findall(r"[-+]?\d*\.\d+|\d+", key)
    
        if value != []:
            data.append(value)
    data = data[5:]
    float_data = [[float(column) for column in row] for row in data]
    float_data = np.array(float_data)
    return float_data

def Xfoil():
    '''
        Расчет профиля крыла(крыла бесконечного удлинения) в Xfoil

        Ввод:   S: float - площадь крыла в плане, м^2
                lambd: float - удлинение крыла
                zeta: float - обратное сужение крыла
                chi_05: float - угол стреловидности по линии середин хорд, рад                
                profile: Profile - объект класса профиля крыла
        Вывод:  ...
    '''

def get_ADX(airfoil, Re, alpha_min, alpha_max, step, ITER):
    '''
        Расчет АДХ профиля (крыла бесконечного удлинения) в Xfoil

        Ввод: airfoil: str - название профиля или имя файла с координатами
              Re: float - число Рейнольдса
              alpha_min: float - минимальное значение угла атаки для расчета, град
              alpha_max: float - максимальное значение угла атаки для расчета, град
              step: float - шаг по углу атаки, град
              ITER: количество итераций численного метода Xfoil
    '''

    try:
        import os
        os.remove('polar1')
        os.remove('polar2')
    except:
        pass

    input_file = open("input_file.in", 'w')
    if ('.dat' in airfoil) or ('.txt' in airfoil):
        input_file.write("load")
        input_file.write('\n')
        input_file.write(airfoil)
        input_file.write('\n')
    else:
        input_file.write(airfoil)
        input_file.write('\n')
    input_file.write("OPER\n")
    input_file.write("visc {0}\n".format(Re))        

    input_file.write("iter {0}\n".format(ITER))
    input_file.write("pacc\n")
    input_file.write("polar1\n")
    input_file.write("polar2\n")
    input_file.write("aseq\n")
    input_file.write("{0}\n".format(alpha_min))
    input_file.write("{0}\n".format(alpha_max))
    input_file.write("{0}\n".format(step))
    input_file.write("quit\n")
    input_file.close()
    subprocess.call("xfoil.exe < input_file.in", shell=True)
    result = {}
    result['alpha'] = read_pack('polar1')[:, 0]
    result['c_ya'] = read_pack('polar1')[:, 1]
    result['c_xa'] = read_pack('polar1')[:, 2]
    result['m_za'] = read_pack('polar1')[:, 4]
    result['x_F'] = - read_pack('polar1')[:, 4] / read_pack('polar1')[:, 1] + 0.25
    os.remove('polar1')
    os.remove('polar2')
    os.remove('input_file.in')
    return result

def get_Cp_x(airfoil, Re, ITER, alpha):
    input_file = open("input_file.in", 'w')
    if ('.dat' in airfoil) or ('.txt' in airfoil):
        input_file.write("load")
        input_file.write('\n')
        input_file.write(airfoil)
        input_file.write('\n')
        input_file.write("mdes")
        input_file.write('\n')
        input_file.write("filt")
        input_file.write('\n')
        input_file.write("exec")
        input_file.write('\n')
        input_file.write('\n')
        input_file.write("pane")
        input_file.write('\n')
    
    else:
        input_file.write(airfoil)
        input_file.write('\n')
    input_file.write("OPER\n")
    input_file.write("visc {0}\n".format(Re))        
    input_file.write("iter {0}\n".format(ITER))
    input_file.write("alfa {0}\n".format(alpha))
    
    input_file.write("cpwr f\n")
    input_file.write("\n\n")
    input_file.write("quit\n")
    input_file.close()
    subprocess.call("xfoil.exe < input_file.in", shell=True)
    result = {}
    result['x'] = read_pack('f')[:, 0]
    result['y'] = read_pack('f')[:, 1]
    result['Cp'] = read_pack('f')[:, 2]
    os.remove('f')
    os.remove('input_file.in')
    return result

def get_airfoil_coords(airfoil_name):
    '''
    Получение координат профиля из Xfoil

    Ввод: airfoil_name: float - название профиля

    Вывод: x: numpy.array - массив координат x профиля
          y: numpy.array - массив координат y профиля
    '''
    
    try:
        import os
        os.remove('polar1')
        os.remove('polar2')
    except:
        pass

    input_file = open("input_file.in", 'w')
    input_file.write(airfoil_name)
    input_file.write('\n')
    input_file.write("psav\n")
    input_file.write("airfoil_coords.dat\n")      
    input_file.write("quit\n")
    input_file.close()
    subprocess.call("xfoil.exe < input_file.in", shell=True)
    
    # чтение созданного Xfoil файла АДХ
    coords = read_pack('airfoil_coords.dat')
    
    os.remove('airfoil_coords.dat')
    os.remove('input_file.in')

    x = coords[:,0]
    y = coords[:,1]

    return x, y

def get_M_cr_Khristianovich(Cp):
    k = 1.4
    h = np.sqrt((k+1) / (k-1))
    def u_def(lam):
        return h*np.sqrt((1- lam**2) / (h**2 - lam**2))
    def L_def(lam):
        return lam*2*np.sqrt(((h-1) ** (h-1)) / ((h+1) ** (h+1))) * np.sqrt(((h+u_def(lam)) ** (h+1)) / (((h - u_def(lam)) ** (h-1)) * ((u_def(lam)+1) ** (2))))
    def a_def(lam):
        return L_def(lam) / lam
    def L_inf_def(Cp):
        return L_def(1) / (np.sqrt(1 - Cp))
    def L_def_for_solve(lam):
        return - L_inf_def(Cp) + lam*2*np.sqrt(((h-1) ** (h-1)) / ((h+1) ** (h+1))) * np.sqrt(((h+u_def(lam)) ** (h+1)) / (((h - u_def(lam)) ** (h-1)) * ((u_def(lam)+1) ** (2))))
    lam_inf_kr = fsolve(L_def_for_solve, 0.1)
    def M_kr_def(lam_inf_kr):
        return lam_inf_kr * (((k+1) / 2) - ((k-1) / 2 * lam_inf_kr**2)) ** (-1/2)
    return M_kr_def(lam_inf_kr)[0]