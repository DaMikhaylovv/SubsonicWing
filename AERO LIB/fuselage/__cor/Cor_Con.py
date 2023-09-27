from math import atan, radians, pi

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM

from .Cor import Cor


class Cor_Con(Cor):
    '''Класс конической кормовой части
    
    Функции:
        c_xa0(M) - расчёт коэффициента силы лобового сопротивления при нулевом угле атаки
        coef_c_x_dn(M, is_active) - расчёт поправки к донному сопротивлению цилиндра
        c_y_alpha(M) - расчёт производной по углу атаки коэффициента нормальной силы, 1/рад
        x_F_alpha(M) - расчёт координаты фокуса по углу атаки (от начала кормовой части), м 
        
        F_surf(x) - расчёт площади поверхности без торцов, м^2
        geometry() - вывод геометрических параметров
        
    '''
    
    def __init__(self, *, D: float, eta_cor: float, lambd_cor: float, D_a: float) -> None:
        # расчет основной геометрии (общей для любой кормовой части)
        super().__init__(D=D, eta_cor=eta_cor, lambd_cor=lambd_cor, D_a=D_a)

        # расчет объема конической кормовой части
        self.W = pi * lambd_cor * D**3 * (1 + eta_cor + eta_cor**2) / 12
 
    # @checkSavedValue
    def c_xa0_p(self, M):
        # Проверяем условие, что угол наклона образующей конуса не более 20 град
        if atan(self.D * (1 - self.eta) / self.L) > radians(20):
            # В данном случае (наклон образующей > 20 градусов) 
            # на всей поверхности кормовой части происходит отрыв потока,
            # поэтому сопротивление собственно кормовой части принимается нулевым,
            # так как вся площадь будет учтена в донном сопротивлении
            return 0
        else:            
            # В данном случае по всей длине кормовой части не происходит отрыв потока
            # вызываем библиотечную функцию для с_xa0 конической кормовой части
            return AeroBDSM.get_Cxcor_Con(M, self.lambd, self.eta)

    def coef_c_x_dn(self, M: float, is_active: bool) -> float:
        # Проверяем условие, что угол наклона образующей конуса не более 20 град
        if atan(self.D * (1 - self.eta) / self.L) > radians(20):            
            # В данном случае (наклон образующей > 20 градусов) 
            # на всей поверхности кормовой части происходит отрыв потока,            
            # поэтому вся площадь миделя учитывается в донном сопротивлении, а не только лишь S_dn
            if is_active == False:  #если сопло не работает учитываем всю площадь  
                return 1
            else:                   #если сопло работает, то учитываем только площадь не занятую соплом
                return (self.S_m - self.S_a) / self.S_m
        else:
            # В данном случае по всей длине кормовой части не происходит отрыв потока
            # поэтому донное сопротивление рассчитывается по фактическому донному срезу кормовой части
            return super().coef_c_x_dn(M, is_active)
        

    def c_y_alpha(self, M: float) -> float:
        return super().c_y_alpha(M)

    def x_F_alpha(self, M: float) -> float:
        return super().x_F_alpha(M)

    # @checkSavedValue
    def F_surf(self, x):
        super().F_surf(x)
        return pi * self.L * (self.D / 2 * (1 + self.eta))