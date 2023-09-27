# Импортирование функции "os.path.dirname"
from os.path import dirname

# Импортирование функции "scipy.interpolate.interp1d"
from scipy.interpolate import interp1d

# Импортирование "pandas"
import pandas as pd

# Импортирование модуля "matplotlib.pyplot"
import matplotlib.pyplot as plt

#  Класс атмосферы "ГОСТ 4401-81 (до 100 км)"
class Atmosphere:
    '''
    Класс атмосферы "ГОСТ 4401-81 (до 100 км)"
    '''
    def __init__(self, datafile_name):
        '''
        Создание полей объекта класса и интерполяция в соответствии с высотой полета

        Функции: get_plot_atm() - функция графического представления параметров атмосферы

        Ввод: datafile_name: str - файл, в котором хранятся значения параметров атмосферы
        Вывод: ...
        '''

        # Наименование колонок данных
        names_of_columns = ['h', 'p', 'rho', 'a', 'g', 'T', 'nu', 'mu']

        # Чтение из .csv файла
        self.DF = pd.read_csv(dirname(__file__) + f'/{datafile_name}', sep = '\t', header = None, names = names_of_columns)

        # Высота
        self.h = self.DF['h']

        # Давление
        self.p_interp = interp1d(self.DF['h'], self.DF['p'])

        # Плотность
        self.rho_interp = interp1d(self.DF['h'], self.DF['rho'])

        # Скорость звука
        self.a_interp = interp1d(self.DF['h'], self.DF['a'])

        # Ускорение свободного падения
        self.g_interp = interp1d(self.DF['h'], self.DF['g'])

        # Температура
        self.T_interp = interp1d(self.DF['h'], self.DF['T'])

        # Коэффициент кинематической вязкости
        self.nu_interp = interp1d(self.DF['h'], self.DF['nu'])

        # Коэффициент динамической вязкости
        self.mu_interp = interp1d(self.DF['h'], self.DF['mu'])

        # Показатель адиабаты воздуха
        self.k = 1.4

        # начальная высота и начальные расчитанные параметры
        self.H = 0.0
        self.saved_params = {}
        self.recalculate(self.H)


    def recalculate(self, H: float) -> dict:
        self.H = H
        self.saved_params = {
            'p': self.p_interp(H),
            'rho': self.rho_interp(H),
            'a': self.a_interp(H),
            'g': self.g_interp(H),
            'T': self.T_interp(H),
            'nu': self.nu_interp(H),
            'mu': self.mu_interp(H)
        }
        return self.saved_params

    def p(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['p']
        self.recalculate(H)
        return self.saved_params['p']

    def rho(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['rho']
        self.recalculate(H)
        return self.saved_params['rho']
    
    def a(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['a']
        self.recalculate(H)
        return self.saved_params['a']
    
    def g(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['g']
        self.recalculate(H)
        return self.saved_params['g']

    def T(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['T']
        self.recalculate(H)
        return self.saved_params['T']

    def nu(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['nu']
        self.recalculate(H)
        return self.saved_params['nu']

    def mu(self, H: float) -> float:
        if self.H == H:
            return self.saved_params['mu']
        self.recalculate(H)
        return self.saved_params['mu']

    def get_plot_atm(self):
        '''
        Графическое представление параметров атмосферы

        Ввод:...
        Вывод:...
        '''

        # Параметры шрифтов графика
        rc = {
            'font.family': 'serif',
            'font.serif': 'Times New Roman',
            'font.variant': 'normal',
            'font.weight': 'normal',
            'font.stretch': 'normal',
            'font.style': 'normal',
            'font.size': 12.0,
            'mathtext.default': 'regular',
            'mathtext.fontset': 'stix'
        }
        plt.rcParams.update(rc)

        # Создание областей графиков
        fig, axs = plt.subplots(2, 4, figsize = (30, 16), sharey = True)

        # Толщина и цвет линий
        w_line = 3
        color_arr = ['blue', 'green', 'orange', 'green', 'lime', 'purple', 'lime', 'grey']

        # Построение графиков
        axs[0][0].plot(self.p(self.h), self.h, c = color_arr[0], linewidth = w_line, alpha = 0.5)
        axs[0][1].plot(self.rho(self.h), self.h, c = color_arr[1], linewidth = w_line, alpha = 0.5)
        axs[0][2].plot(self.T(self.h), self.h, c = color_arr[2], linewidth = w_line, alpha = 0.5)
        axs[0][3].plot(self.a(self.h), self.h, c = color_arr[3], linewidth = w_line, alpha = 0.5)
        axs[1][0].plot(self.g(self.h), self.h, c = color_arr[4], linewidth = w_line, alpha = 0.5)
        axs[1][1].plot(self.nu(self.h), self.h, c = color_arr[5], linewidth = w_line, alpha = 0.5)
        axs[1][2].plot(self.mu(self.h), self.h, c = color_arr[6], linewidth = w_line, alpha = 0.5)
        axs[1][3].plot(self.h, self.h, c = color_arr[7], linewidth = w_line, alpha = 0.5)

        # Линии сетки
        axs[0][0].grid()
        axs[0][1].grid()
        axs[0][2].grid()
        axs[0][3].grid()
        axs[1][0].grid()
        axs[1][1].grid()
        axs[1][2].grid()
        axs[1][3].grid()

        # Подписи
        axs[0][0].set_title(r'$Давление$')
        axs[0][0].set_xlabel(r"$p$, $Па$")
        axs[0][0].set_ylabel(r"$h$, $м$")
        axs[0][1].set_title(r'$Плотность$')
        axs[0][1].set_xlabel(r"$\rho$, $кг/м^3$")
        axs[0][2].set_title(r'$Температура$')
        axs[0][2].set_xlabel(r"$T$, $K$")
        axs[0][3].set_title(r'$Скорость$ $звука$')
        axs[0][3].set_xlabel(r"$a$, $м/с$")
        axs[1][0].set_title(r'$Ускорение$ $свободного$ $падения$')
        axs[1][0].set_xlabel(r"$g$, $м/с^2$")
        axs[1][0].set_ylabel(r"$h$, $м$")
        axs[1][1].set_title(r'$Коэффициент$ $кинематической$ $вязкости$')
        axs[1][1].set_xlabel(r"$\nu$, $м^2/с$")
        axs[1][2].set_title(r'$Коэффициент$ $динамической$ $вязкости$')
        axs[1][2].set_xlabel(r"$\mu$, $Па∙с$")
        axs[1][2].set_title(r'$Высота$')
        axs[1][3].set_xlabel(r"$h$, $м$")
        
        # Сохранение графика
        plt.savefig('ГОСТ 4401-81 (до 100 км).png')
        plt.show()

# Переменная имени файла с данными атмосферы
__atmosphere_filename = 'atmosphere_GOST_data.csv'

# Объект класса атмосферы
atmo = Atmosphere(__atmosphere_filename)
