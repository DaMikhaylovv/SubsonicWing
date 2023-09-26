# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Import necessary modules
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import re
import os.path
import numpy as np
from scipy.optimize import fsolve
import subprocess
from scipy.interpolate import interp1d

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           Core Functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Session:
    VALUE = {'Re':[50e3, 8e6], 'ITER':[20, 250], 'INIT_TYPE': ['alpha', 'Cl'], 'alpha': [-25, 25], 'Cl': [-2, 2]}
    AIRFOILS = ['NACA2208', 'NACA2211', 'NACA2214', 'NACA2217', 
                'NACA2220', 'NACA23008', 'NACA23011', 'NACA23014', 
                'NACA23017', 'NACA23020', 'NACA0012']
    METHOD_M_KR = ['KHRISTIANOVICH', 'BURAGO', 'KARMAN']

    def get_available_airfoils(self):
        """ 
        Возвращает список доступных профилей (уточняется)

            Выходные параметры
            ==================
        
            AIRFOIL (list): список доступных профилей
        """
        return self.AIRFOILS
    def get_available_method_M_kr(self):
        """ 
        Возвращает список доступных методов расчета критического числа Маха

            Выходные параметры
            ==================
            METHOD_M_KR (list): список доступных профилей
        """
        return self.METHOD_M_KR


    @classmethod
    def check_value(cls, value, key):
        # if key == 'increment':
        #     return type(value) in (float, int)
        if key == 'AIRFOIL':
            if type(value) is str:
                return value in cls.AIRFOILS
            else:
                return False
        elif key == 'AIRFOIL_FILE':
            return type(value) is str
        elif key == 'Re':
            if type(value) in (float, int):
                return (cls.VALUE['Re'][0] <= value <= cls.VALUE['Re'][1])
            else:
                return False
        elif key == 'ITER':
            if type(value) in (float, int):
                return cls.VALUE['ITER'][0] <= value <= cls.VALUE['ITER'][1]
            else:
                return False
        elif key == 'INIT_TYPE':
            return value in cls.VALUE['INIT_TYPE']
        elif key == 'alpha':
            if type(value) in (float, int):
                return cls.VALUE['alpha'][0] <= value <= cls.VALUE['alpha'][1]
            else:
                return False
        elif key == 'Cl':
            if type(value) in (float, int):
                return cls.VALUE['Cl'][0] <= value <= cls.VALUE['Cl'][1]
            else:
                return False
        elif key == 'method':
            if type(value) is str:
                return value in cls.METHOD_M_KR
            else:
                return False
    @staticmethod
    def check_file_loc(file):
        return os.path.exists(file)

    @staticmethod
    def read_pack(filename):
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

    @staticmethod
    def Khristianovich(Cp):
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
    
    @staticmethod
    def Burago_by_Fabricant(Cp):
        k = 1.4
        def p_solve(M_kr):
            return - Cp + 1 - (1 / M_kr**2) * ((2 / (k+1)) + ((k-1) / (k+1))*M_kr**2)**(k / (k-1))
        M_kr = fsolve(p_solve, 0.5)
        return M_kr[0]
    
    @staticmethod
    def Karman_2(Cp):
        k=1.4
        def for_solve(M_kr):
            return -Cp + 2 * np.sqrt(1-M_kr**2) * (1 - (((2 + (k-1)*M_kr**2) / (k+1)) ** (k / (k-1)))) / (k*M_kr**2)
        res = fsolve(for_solve, 0.5)[0]
        return res

    @staticmethod
    def read_pack(filename):
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

    # def __init__(self, airfoil, Re, ITER, INIT_TYPE, value):
    #     if ('.dat' in airfoil) or ('.txt' in airfoil):
    #         if self.check_value(airfoil, 'AIRFOIL_FILE'):
    #             if self.check_file_loc(airfoil):
    #                 self.airfoil = airfoil
    #             else:
    #                 raise NameError("Файл отсутствует")
    #         else:
    #             raise ValueError("Неверное имя файла с координатами профиля")
    #     else:
    #         if self.check_value(airfoil, 'AIRFOIL'):
    #             self.airfoil = airfoil
    #         else:
    #             raise ValueError("Недоступный профиль или неверный формат")

    #     if self.check_value(Re, 'Re'):
    #         self.Re = Re
    #     else:
    #         raise ValueError("Недопустимое значение числа Рейнольдса")
    #     if self.check_value(ITER, 'ITER'):
    #         self.ITER = ITER
    #     else:
    #         raise ValueError("Недопустимое значение числа итераций")
    #     if self.check_value(INIT_TYPE, 'INIT_TYPE'):
    #         self.INIT_TYPE = INIT_TYPE
    #     else:
    #         raise ValueError("Неверно задана переменная INIT_TYPE")
    #     if INIT_TYPE == 'alpha':
    #         if self.check_value(value, 'alpha'):
    #             self.alpha = value
    #         else:
    #             raise ValueError("Недопустимое значение угла атаки")
    #     else:
    #         if self.check_value(value, 'Cl'):
    #             self.Cl = value
    #         else:
    #             raise ValueError("Недопустимое значение коэффициента подъемной силы")

    def get_Cp_x(self, airfoil, Re, ITER, INIT_TYPE, value):
        """
        Аргументы
        =========

            airfoil: 'str'
                профиль крыла;
            Re: 'float'
                число Рейнольдса [50000 ... 2000000];
            ITER: 'int':
                число итераций [20...250];
            INIT_TYPE: 'str':
                переменная инициализирующая расчет по углу атаки ('alpha') или коэффициенту
                подъемной силы c_ya ('Cl');
            value: 'float':
                значение угла атаки alpha (если INIT_TYPE == 'alpha') или коэффициента
                подъемной силы c_ya (есди INIT_TYPE == 'Cl')

        Выходные параметры
        ==================

        result: 'dict'
            словарь следующего вида
            {['x']: list (координата x профиля),
            ['y']: list (координата y профиля),
            ['Cp']: list (значение коэффициента давления Cp)
            }
        """

        if ('.dat' in airfoil) or ('.txt' in airfoil):
            if self.check_value(airfoil, 'AIRFOIL_FILE'):
                if self.check_file_loc(airfoil):
                    self.airfoil = airfoil
                else:
                    raise NameError("Файл отсутствует")
            else:
                raise ValueError("Неверное имя файла с координатами профиля")
        else:
            if self.check_value(airfoil, 'AIRFOIL'):
                self.airfoil = airfoil
            else:
                raise ValueError("Недоступный профиль или неверный формат")

        if self.check_value(Re, 'Re'):
            self.Re = Re
        else:
            raise ValueError("Недопустимое значение числа Рейнольдса")
        if self.check_value(ITER, 'ITER'):
            self.ITER = ITER
        else:
            raise ValueError("Недопустимое значение числа итераций")
        if self.check_value(INIT_TYPE, 'INIT_TYPE'):
            self.INIT_TYPE = INIT_TYPE
        else:
            raise ValueError("Неверно задана переменная INIT_TYPE")
        if INIT_TYPE == 'alpha':
            if self.check_value(value, 'alpha'):
                self.alpha = value
            else:
                raise ValueError("Недопустимое значение угла атаки")
        else:
            if self.check_value(value, 'Cl'):
                self.Cl = value
            else:
                raise ValueError("Недопустимое значение коэффициента подъемной силы")
        input_file = open("input_file.in", 'w')
        if ('.dat' in self.airfoil) or ('.txt' in self.airfoil):
            input_file.write("load")
            input_file.write('\n')
            input_file.write(self.airfoil)
            input_file.write('\n')
            # input_file.write("airfoil")
            # input_file.write('\n')
            input_file.write("mdes")
            input_file.write('\n')
            input_file.write("filt")
            input_file.write('\n')
            input_file.write("exec")
            input_file.write('\n')
            input_file.write('\n')
            input_file.write("pane")
            input_file.write('\n')
            # input_file.write("{0}\n".format(400))
            # input_file.write('\n')

        else:
            input_file.write(self.airfoil)
            input_file.write('\n')
        input_file.write("OPER\n")
        input_file.write("visc {0}\n".format(self.Re))        
        input_file.write("iter {0}\n".format(self.ITER))
        
        if self.INIT_TYPE == 'alpha':
            input_file.write("alfa {0}\n".format(self.alpha))
        elif self.INIT_TYPE == 'Cl':
            input_file.write("Cl {0}\n".format(self.Cl))
        
        input_file.write("cpwr f\n")
        input_file.write("\n\n")
        input_file.write("quit\n")
        input_file.close()
        subprocess.call("xfoil.exe < input_file.in", shell=True)
        result = {}
        result['x'] = self.read_pack('f')[:, 0]
        result['y'] = self.read_pack('f')[:, 1]
        result['Cp'] = self.read_pack('f')[:, 2]
        return result        

    def get_Cp_min(self, Cp):
        """
        Функция для расчета минимального значения коэффициента давления

        Аргементы
        =========

            Cp: 'list'
                Список значений коэффициентов давлений

        Выходные параметры
        ==================

            Cp_min: 'float'
                Минимальное значение коэффициента давления Cp
        """

        return np.min(Cp)

    def get_M_kr(self, Cp_min, method):
        """
        Функция для расчета критического числа Маха

        Аргументы
        =========

            Cp_min: 'float'
                Минимальное значение коэффициента давления
            method: 'str'
                Метод расчета критического числа Маха

        Выходные параметры
        ==================

            M_kr: 'float'
                Критическое число Маха
        """

        if self.check_value(method, 'method'):
            if method == 'KHRISTIANOVICH':
                return self.Khristianovich(Cp_min)
            elif method == 'KHRISTIANOVICH':
                return self.Burago_by_Fabricant(Cp_min)
            elif method == 'BURAGO':
                return self.Karman_2(Cp_min)
        else:
            raise NameError("Неверно задан метод расчета M_кр")

    def get_ADX(self, airfoil, Re, forADX, value_min, value_max, inc, ITER):
        """
        Аргументы
        =========

            airfoil: 'str'
                профиль крыла;
            Re: 'float'
                число Рейнольдса [50000 ... 2000000];
            forADX: 'str':
                переменная, инициализирующая расчет АДХ профиля крыла по по углу атаки ('alpha') 
                или коэффициенту подъемной силы c_ya ('Cl');
            value_min: 'float':
                начальное значение угла атаки alpha (если INIT_TYPE == 'alpha') или коэффициента
                подъемной силы c_ya (есди INIT_TYPE == 'Cl');
            value_max: 'float':
                конечное значение угла атаки alpha (если INIT_TYPE == 'alpha') или коэффициента
                подъемной силы c_ya (есди INIT_TYPE == 'Cl');
            inc: 'float':
                величина шага;
            ITER: 'int':
                число итераций [20...180];

        Выходные параметры
        ==================

        result: 'dict'
            словарь следующего вида
            {['alpha']: list (углы атаки),
            ['c_ya']: list (коэффициенты подъемной силы c_ya),
            ['c_xa']: list (коэффициенты силы сопротивления c_xa)
            }
        """
        try:
            import os
            os.remove('polar1')
            os.remove('polar2')
        except:
            pass
        if forADX == 'alpha' or forADX == 'Cl':
            if forADX == 'alpha':
                if self.check_value(value_min, 'alpha') and self.check_value(value_max, 'alpha'):
                    self.alpha_min = value_min
                    self.alpha_max = value_max
                    self.inc = inc
                else:
                    raise ValueError("Недопустимое значение угла атаки")
            elif forADX == 'Cl':
                if self.check_value(value_min, 'Cl') and self.check_value(value_max, 'Cl'):
                    self.Cl_min = value_min
                    self.Cl_max = value_max
                    self.inc = inc
                else:
                    raise ValueError("Недопустимое значение коэффициента подъемной силы")
        else:
            raise ValueError("Недопустимое задание исходных данных для расчета")
        
        # elif forADX == 'Cl':
        #     if self.check_value(value, 'Cl'):
        #         self.Cl_ADX = value
        #     else:
        #         raise ValueError("Недопустимое значение коэффициента подъемной силы")
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
        if forADX == 'alpha':
            input_file.write("aseq\n")
            input_file.write("{0}\n".format(self.alpha_min))
            input_file.write("{0}\n".format(self.alpha_max))
            input_file.write("{0}\n".format(self.inc))
            # input_file.write("alfa {0}\n".format(self.alpha_ADX))
        elif forADX == 'Cl':
            input_file.write("cseq\n")
            input_file.write("{0}\n".format(self.Cl_min))
            input_file.write("{0}\n".format(self.Cl_max))
            input_file.write("{0}\n".format(self.inc))
            # input_file.write("Cl {0}\n".format(self.Cl_ADX))
        input_file.write("quit\n")
        input_file.close()
        subprocess.call("xfoil.exe < input_file.in", shell=True)
        result = {}
        result['alpha'] = self.read_pack('polar1')[:, 0]
        result['c_ya'] = self.read_pack('polar1')[:, 1]
        result['c_xa'] = self.read_pack('polar1')[:, 2]
        os.remove('polar1')
        os.remove('polar2')
        return result

class Experiment:
    AIRFOILS = ['NACA2208', 'NACA2211', 'NACA2214', 'NACA2217', 
                'NACA2220', 'NACA23008', 'NACA23011', 'NACA23014', 
                'NACA23017', 'NACA23020', 'NACA0012']
    
    @classmethod
    def check_value(cls, value, key):
        if key == 'AIRFOIL':
            if type(value) is str:
                return value in cls.AIRFOILS
            else:
                return False

    def get_M_kr_Cp_values(self):
        """
        Функция, возвращающая результаты значений Cp_min и M_kr [Атлас, 1940 г.]

        Выходные параметры
        ==================

        result: 'dict'
            словарь следующего вида
            {['Cp_min']: list (значения минимальных коэффициентов давления Cp_min),
            ['M_kr']: list (значения критических чисел Маха)
            }
        """
        Cp_min = np.array([0, -0.10, -0.20, -0.30, -0.40, -0.50, -0.60, -0.70, -0.80,
                        -0.90, -1.00, -1.25, -1.50, -1.75, -2.00, -2.25, -2.50, -2.60,
                        -2.70, -2.80, -2.90, -3.00])
        M_kr = np.array([1, 0.889, 0.811, 0.750, 0.704, 0.667, 0.636, 0.610,
                         0.588, 0.568, 0.550, 0.511, 0.479, 0.454, 0.433, 0.413,
                         0.395, 0.388, 0.382, 0.376, 0.371, 0.366])
        data = {}
        data['Cp_min'] = Cp_min
        data['M_kr'] = M_kr
        return data

    def get_Cp_experiment(self, airfoil):
        """
        Функция, возвращающая зависимости Cp(x), полученные экспериментально [Атлас, 1940 г.]

        Аргументы
        =========

            airfoil: 'str'
                наименование профиля крыла

        Выходные параметры
        ==================

        Cp_x_up: numpy.ndarray:
            структура массива как в атласе 1940 года (стр. 186)
        Cp_x_up: numpy.ndarray:
            структура массива как в атласе 1940 года (стр. 186)
            
        """
        NACA_0012_up = np.array([[1, 0.7, 0.210, -0.180, -0.570, -1.020, -1.480, -2.420, -3.340],
                         [2.5, 0, -0.210, -0.400, -0.660, -0.955, -1.270, -1.990, -2.850],
                         [5.0, -0.2, -0.35, -0.515, -0.680, -0.885, -1.090, -1.620, -2.150],
                         [10, -0.3, -0.4, -0.525, -0.640, -0.790, -0.930, -1.235, -1.550],
                         [15, -0.3, -0.390, -0.50, -0.60, -0.715, -0.810, -1.040, -1.250],
                         [20, -0.290, -0.370, -0.465, -0.550, -0.645, -0.715, -0.890, -1.060],
                         [30, -0.260, -0.320, -0.390, -0.450, -0.515, -0.565, -0.680, -0.805],
                         [40, -0.220, -0.270, -0.310, -0.355, -0.405, -0.440, -0.525, -0.610],
                         [50, -0.175, -0.210, -0.235, -0.270, -0.300, -0.325, -0.390, -0.440],
                         [60, -0.130, -0.140, -0.160, -0.190, -0.205, -0.220, -0.260, -0.285],
                         [70, -0.080, -0.090, -0.100, -0.120, -0.115, -0.120, -0.145, -0.140],
                         [80, -0.040, -0.050, -0.040, -0.045, -0.035, -0.040, -0.040, -0.020],
                         [90, -0.010, -0.010, 0.015, 0.030, 0.040, 0.050, 0.060, 0.075],
                         [95, 0.005, 0.020, 0.040, 0.065, 0.070, 0.080, 0.095, 0.105]])
        NACA_0012_dn = np.array([[1, 0, 0.210, 0.440, 0.580, 0.680, 0.820, 0.960, 1.00], 
                         [2.5, -0.380, -0.210, 0.020, 0.230, 0.325, 0.510, 0.760, 0.950],
                         [5, -0.5, -0.35, -0.190, -0.015, 0.110, 0.260, 0.530, 0.730],
                         [10, -0.520, -0.400, -0.300, -0.175, -0.055, 0.060, 0.280, 0.450], 
                         [15, -0.485, -0.390, -0.310, -0.215, -0.120, -0.030, 0.155, 0.310],
                         [20, -0.445, -0.370, -0.300, -0.215, -0.140, -0.075, 0.090, 0.225],
                         [30, -0.370, -0.320, -0.270, -0.200, -0.140, -0.100, 0.020, 0.135],
                         [40, -0.3, -0.270, -0.210, -0.175, -0.125, -0.080, 0.005, 0.100],
                         [50, -0.240, -0.210, -0.155, -0.140, -0.100, -0.060, 0.015, 0.080],
                         [60, -0.180, -0.140, -0.105, -0.090, -0.065, -0.030, 0.025, 0.070],
                         [70, -0.120, -0.090, -0.060, -0.040, -0.025, -0.005, 0.035, 0.060],
                         [80, -0.060, -0.050, -0.020, -0.010, 0.000, 0.015, 0.040, 0.055],
                         [90, -0.010, -0.010, 0.010, 0.015, 0.025, 0.025, 0.035, 0.040],
                         [95, 0.030, 0.020, 0.025, 0.030, 0.040, 0.040, 0.060, 0.030]])

        NACA_23008_up = np.array([[1, 0.790, 0.660, 0.400, 0.170, -0.110, -0.400, -1.240, -2.400],
                         [2.5, 0.490, 0.290, 0.060, -0.170, -0.410, -0.690, -1.350, -2.160],
                         [5.0, 0.120, -0.070, -0.230, -0.420, -0.635, -0.860, -1.360, -1.910],
                         [10,  -0.240, -0.365, -0.520, -0.650, -0.780, -0.930, -1.210, -1.550],
                         [15,  -0.270, -0.370, -0.490, -0.595, -0.710, -0.800, -1.020, -1.240],
                         [20,  -0.220, -0.310, -0.400, -0.490, -0.570, -0.660, -0.820, -0.990],
                         [30,  -0.175, -0.240, -0.285, -0.330, -0.400, -0.430, -0.570, -0.680],
                         [40,  -0.145, -0.185, -0.215, -0.270, -0.300, -0.340, -0.425, -0.510],
                         [50,  -0.120, -0.140, -0.160, -0.200, -0.225, -0.260, -0.320, -0.380],
                         [60,  -0.090, -0.100, -0.120, -0.140, -0.165, -0.190, -0.225, -0.270],
                         [70,  -0.050, -0.065, -0.080, -0.090, -0.100, -0.115, -0.145, -0.170],
                         [80,  -0.010, -0.030, -0.035, -0.035, -0.040, -0.050, -0.065, -0.070],
                         [90,  0.040, 0.030, 0.030, 0.035, 0.020, 0.020, 0.010, 0.020],
                         [95,  0.080, 0.080, 0.070, 0.080, 0.070, 0.070, 0.060, 0.060]])
        NACA_23008_dn = np.array([[1, -1.150, -0.820, -0.500, -0.180, 0.100, 0.400, 0.800, 0.990], 
                         [2.5,  -1.020, -0.710, -0.430, -0.190, 0.040, 0.240, 0.600, 0.800],
                         [5, -0.540, -0.385, -0.200, -0.050, 0.100, 0.230, 0.490, 0.680],
                         [10, -0.264, -0.160, -0.050, 0.050, 0.150, 0.235, 0.410, 0.560], 
                         [15, -0.230, -0.150, -0.060, 0.020, 0.100, 0.180, 0.325, 0.450],
                         [20, -0.290, -0.205, -0.130, -0.045, 0.020, 0.095, 0.235, 0.360],
                         [30, -0.250, -0.190, -0.130, -0.075, 0.000, 0.040, 0.150, 0.250],
                         [40, -0.190, -0.150, -0.100, -0.060, 0.000, 0.030, 0.120, 0.200],
                         [50, -0.140, -0.110, -0.070, -0.040, 0.000, 0.030, 0.100, 0.160],
                         [60, -0.095, -0.070, -0.040, -0.020, 0.005, 0.035, 0.090, 0.135],
                         [70, -0.060, -0.040, -0.020, 0.000, 0.020, 0.040, 0.080, 0.120],
                         [80, -0.030, 0.000, 0.000, 0.020, 0.035, 0.050, 0.080, 0.100],
                         [90, 0.040, 0.000, 0.000, 0.020, 0.035, 0.050, 0.080, 0.100],
                         [95, 0.080, 0.085, 0.080, 0.090, 0.095, 0.100, 0.110, 0.100]])

        if self.check_value(airfoil, 'AIRFOIL'):
            if airfoil == 'NACA0012':
                return NACA_0012_up, NACA_0012_dn
            
            elif airfoil == 'NACA23008':
                return NACA_23008_up, NACA_23008_dn
        else:
            raise ValueError("Неверно указан профиль")

