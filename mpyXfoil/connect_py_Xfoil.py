import os.path
import subprocess
import numpy as np
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
    return result