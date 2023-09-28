from abc import ABC
from math import *

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM


class Cor(ABC):
    '''Интерфейс классов кормовой части
    
    Функции:
        c_xa0(M) - расчёт коэффициента силы лобового сопротивления при нулевом угле атаки
        coef_c_x_dn(M, is_active) - расчёт поправки к донному сопротивлению цилиндра
        c_y_alpha(M) - расчёт производной по углу атаки коэффициента нормальной силы, 1/рад
        x_F_alpha(M) - расчёт координаты фокуса по углу атаки (от начала кормовой части), м 
        
        F_surf(x) - расчёт площади поверхности без торцов, м^2
        geometry() - вывод геометрических параметров
        
    '''
    def __init__(self, *, D: float, eta_cor: float, lambd_cor: float, D_a: float) -> None:
        self.D = D                                  # Диаметр миделя
        self.eta = eta_cor                          # Сужение
        self.lambd = lambd_cor                      # Удлинение
        self.D_a = D_a                              # Диаметр выходного сечения сопла

        self.D_dn = self.eta * self.D               # Диаметр донного среза
        self.L = self.D * self.lambd                # Длина        
        self.S_m = pi * self.D**2 / 4               # Площадь миделя
        self.S_dn = pi * self.D_dn**2 / 4           # Площадь донного среза
        self.S_a = pi * self.D_a**2 / 4             # Площадь выходного сечения сопла
        self.W = None                               # Объем        
        
    def c_xa0_p(self, M: float) -> float:
        '''
        Расчёт коэффициента силы сопротивления давления при нулевом угле атаки

        Ввод:   M: float - число Маха

        Вывод:  c_xa0 - коэффициент силы лобового сопротивления при нулевом угле атаки
        '''
        ... 

    # @checkSavedValue
    def coef_c_x_dn(self, M: float, is_active: bool) -> float:
        '''
        Расчёт поправки к донному сопротивлению цилиндра        

        Ввод:   M: float - число Маха
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  coef_c_x_dn - поправка к донному сопротивлению цилиндра, учитывающая сужение кормовой части и сопло      
        '''
        # Вызываем библиотечную функцию для определения поправочного коэффициента, учитывающего сужение кормовой части        
        Sigma_eta = AeroBDSM.get_Sigma_eta(M, self.lambd, self.eta)
        
        # Вводим поправку на то, что площадь донного среза кормовой части отличается от площади цилиндра
        if is_active == False:  #если сопло не работает учитываем всю площадь донного среза  
            return Sigma_eta * self.S_dn / self.S_m
        else:                   #если сопло работает, то учитываем только площадь дна не занятую соплом
            return Sigma_eta * (self.S_dn - self.S_a) / self.S_m
            
    # @checkSavedValue
    def c_y_alpha(self, M: float) -> float:
        '''
        Расчёт производной коэффициента нормальной силы по углу атаки
        
        Ввод:   M: float - число Маха
                        
        Вывод:  c_y_alpha: float - производная коэффициента нормальной силы по углу атаки, 1/рад
        '''
        return -0.2 * 2 * (1 - self.eta**2)
    
    # @checkSavedValue
    def x_F_alpha(self, M: float) -> float:
        '''
        Расчёт координаты фокуса по углу атаки

        Ввод:   M: float - число Маха
                        
        Вывод:  x_F_alpha: float - координата фокуса по углу атаки (от начала кормовой части), м
        '''
        
        return self.L - (self.S_m * self.L - self.W) / (self.S_m - pi * self.D_dn**2 / 4)
        

    def F_surf(self, x):
        '''
        Расчёт площади поверхности
        
        Ввод:   x: float - координата границы области расчёта площади поверхности (от начала кормовой части), м
        
        Вывод:  F_surf: float - площадь поверхности (без торцов), м^2
        '''
        # Проверяем правильность входных данных
        assert all([x >= 0, x <= self.L]), \
            ValueError(f'Используемая координата выходит за пределы кормовой части, значение координаты: {x}')

    def geometry(self):
        '''
        Функция вывода геометрии.

        Вывод: Dict[str, float], где ключи -- названия геометрического параметра
        '''
        attrlist = ['D', 'eta', 'lambd', 'D_a', 'D_dn', 'L', 'S_m', 'S_dn', 'S_a', 'W']
        
        return {attrname:self.__getattribute__(attrname) for attrname in attrlist}