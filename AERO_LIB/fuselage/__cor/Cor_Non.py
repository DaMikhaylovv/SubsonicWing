from math import pi
from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM

from .Cor import Cor


class Cor_Non(Cor):
    '''Класс отсутствующей кормовой части
    
    Функции:
        c_xa0(M) - расчёт коэффициента силы лобового сопротивления при нулевом угле атаки
        coef_c_x_dn(M, is_active) - расчёт поправки к донному сопротивлению цилиндра
        c_y_alpha(M) - расчёт производной по углу атаки коэффициента нормальной силы, 1/рад
        x_F_alpha(M) - расчёт координаты фокуса по углу атаки (от начала кормовой части), м 
        
        F_surf(x) - расчёт площади поверхности без торцов, м^2
        geometry() - вывод геометрических параметров
        
    '''


    def __init__(self, *, D: float, D_a: float) -> None:
        # Расчет основной геометрии (общей для любой кормовой части)
        super().__init__(D=D, eta_cor=1, lambd_cor=0, D_a=D_a)
        
        # Расчет объема кормовой части
        self.W = 0   

    def c_xa0_p(self, M):
        return 0
    
    def coef_c_x_dn(self, M: float, is_active: bool) -> float:
        if is_active == False:  #если сопло не работает учитываем всю площадь дна
            return 1
        else:                   #если сопло работает, то учитываем только площадь не занятую соплом
            return (self.S_m - self.S_a) / self.S_m
    
    def c_y_alpha(self, M: float) -> float:
        return 0

    def x_F_alpha(self, M: float) -> float:
        return 0

    def F_surf(self, x):
        return 0
