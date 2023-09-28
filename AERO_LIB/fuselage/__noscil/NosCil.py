from math import pi
from abc import ABC
from scipy.interpolate import interp1d

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM


class NosCil(ABC):
    '''
    Интерфейс классов носовой части   
    '''
    def __init__(self, *, lambd_nos: float, D: float, lambd_cil: float, r_sph: float = 0) -> None:
        '''
        Создание геометрии класса:
            D - диаметр миделя, м
            lambd_cil - удлинение цилиндрической части
            lambd_nos - удлинение носовой части
        '''
      
        # Инициализация геометрии
        self.D = D                              # Диаметр миделя
        self.lambd_nos = lambd_nos              # Удлинение носовой части
        self.lambd_cil = lambd_cil              # Удлинение цилиндрической части
        self.r_sph = r_sph                      # Радиус сферического затупления (для сфера-конус и сфера-оживало)        
        
        self.L_nos = self.D * self.lambd_nos    # Длина носовой части
        self.L_cil = self.D * self.lambd_cil    # Длина цилиндрической части
        self.L = self.L_cil + self.L_nos        # Общая длина носовой+цилиндрической части
        self.rr_sph = 2 * self.r_sph / self.D   # Относительный радиус сферического затупления
        self.S_m = pi * D**2 / 4               # Площадь миделя
        
        self.S_bok_nos = None                   # Площадь боковой проекции носовой части
        self.W_nos = None                       # Объем носовой части        
        self.x_cs_nos = None                    # Координата геометрического центра тяжести сечения носовой части


    def c_xa0_p(self, M: float, H: float) -> float:
        '''
        Расчёт коэффициента силы лобового сопротивления давления носовой части при нулевом угле атаки

        Ввод:   M: float - число Маха
                H: float - высота полета, м

        Вывод:  c_xa0 - коэффициент силы лобового сопротивления давления носовой части при нулевом угле атаки
        '''        
        ...

    def Delta_c_x(self, M: float, alpha: float) -> float:
        '''
        Расчёт коэффициента дополнительной продольной силы, возникающей при ненулевых углах атаки,
        в результате перераспределения давления на носовой части

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                        
        Вывод:  Delta_c_x: float - коэффициент дополнительной продольной силы
        '''
        ...

    def c_y_alpha(self, M: float) -> float:
        '''
        Расчёт производной коэффициента нормальной силы по углу атаки
        
        Ввод:   M: float - число Маха
                        
        Вывод:  c_y_alpha: float - производная коэффициента нормальной силы по углу атаки, 1/рад
        '''
        ...
    
    # @checkSavedValue
    def x_F_alpha(self, M: float) -> float:
        '''
        Расчёт координаты фокуса по углу атаки (от носа)

        Ввод:   M: float - число Маха [-]
                        
        Вывод:  x_F_alpha: float - координата фокуса по углу атаки (от носа), м
        '''
        
        # Вызываем библиотечную функцию для определения смещения фокуса с учётом М относительно фокуса по теории тонких тел
        Delta_x_F = AeroBDSM.get_xi_M(M, self.lambd_nos, self.lambd_cil) * self.L_nos

        return self.L_nos - self.W_nos / self.S_m + Delta_x_F

  
    def lambd_nos_eqvCon(self, M: float) -> float:
        '''
        Расчёт удлинения эквивалентного конуса
        
        Ввод:   M: float - число Маха
                        
        Вывод:  lambd_nos_eqvCon: float - удлинение эквивалентной конической носовой части
        '''
        ...

    # @checkSavedValue
    def kM(self, M: float) -> float:
        '''
        Расчёт коэффициента торможения потока, вызванного обтеканием носовой части.
        
        Ввод:   M: float - число Маха [-]
                        
        Вывод:  kM: float - коэффициент торможения потока, вызванного обтеканием носовой части
        '''
        return AeroBDSM.get_kM_Con(M, self.lambd_nos_eqvCon(M)).Value
    
    
    def F_surf(self, x: float) -> float:
        '''
        Расчёт площади поверхности без торцов
        
        Ввод:   x: float - координата границы области расчёта площади поверхности, м
        
        Вывод:  F_surf: float - площадь боковой поверхности без торцов, м^2
        '''
        # Проверяем правильность входных данных
        assert all([x >= 0, x <= self.L]), \
            ValueError(f'Используемая координата выходит за пределы носовой и цилиндрической части фюзеляжа, значение координаты: {x}')
        
    def _x_cs_nos(self, r: callable, step = 1000):
        '''
        Расчет координаты геометрического центра тяжести площади боковой проекции носовой части ??а аналитически для каждой носовой нельзя найти?

        Ввод:   r: callable - функция радиуса поверхности (образующей) от координаты от носа, м
                step: int - (опционально) количество шагов разбиения

        Вывод: x_cs_nos: float - координата тяжести площади боковой проекции носовой части (от носа), м
        '''
        # Расчёт диаметра носовой части в зависимости от координаты Х
        D = lambda x: r(x) * 2

        # Расчёт элемента суммы знаменателя функции геометрического центра тяжести
        denom = lambda i, dx = self.L_nos / step: (D(dx * i) + D(dx * i + dx)) * dx / 2

        # Расчёт элемента суммы числителя функции геометрического центра тяжести
        nom = lambda i, dx = self.L_nos / step: denom(i) * (dx * i + dx / 3 * (2 * D(dx * i + dx) + D(dx * i)) / (D(dx * i) + D(dx * i + dx)))

        # Расчёт геометрического центра тяжести
        self.x_cs_nos = sum(map(nom, range(step))) / sum(map(denom, range(step)))
        return self.x_cs_nos

    def geometry(self):
        '''
        Функция вывода геометрии.

        Вывод: Dict[str, float], где ключи -- названия геометричемского параметра
        '''
        attrlist = [
            'D', 'lambd_nos', 'lambd_cil', 'r_sph', 'S_bok_nos', 'W_nos', 'theta_con', 'x_cs_nos', 'L_cil', 'L_nos', 'S_m', 'L', 'rr_sph']
        if hasattr(self, 'L_main'):
            attrlist.append('L_main')
        return {attrname:self.__getattribute__(attrname) for attrname in attrlist}
