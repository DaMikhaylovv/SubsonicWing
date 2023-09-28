from math import atan, radians, pi, sqrt
import scipy.integrate as integrate

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM


from .Cor import Cor


class Cor_Par(Cor):
    '''Класс параболической кормовой части
    
    Функции:
        c_xa0(M) - расчёт коэффициента силы лобового сопротивления при нулевом угле атаки
        coef_c_x_dn(M, is_active) - расчёт поправки к донному сопротивлению цилиндра
        c_y_alpha(M) - расчёт производной по углу атаки коэффициента нормальной силы, 1/рад
        x_F_alpha(M) - расчёт координаты фокуса по углу атаки (от начала кормовой части), м 
        
        F_surf(x) - расчёт площади поверхности без торцов, м^2
        geometry() - вывод геометрических параметров
        
    '''


    def __init__(self, *, D: float, eta_cor: float, lambd_cor: float, D_a: float) -> None:
        # Расчет основной геометрии (общей для любой кормовой части)
        super().__init__(D=D, eta_cor=eta_cor, lambd_cor=lambd_cor, D_a=D_a)

        # Расчет объема кормовой части
        self.W = integrate.quad(lambda x: (D / 2 * (1 - self.eta) * (x / self.L)**2)**2, 0, self.L)[0] * pi

    def c_xa0_p(self, M):
        # Уравнение образующей кормовой части
        r = lambda x: self.D / 2 * (1 - self.eta) * (x / self.L)**2

        # Рассчитываем функцию угла наклона касательной к образующей  ??а нельзя это аналитически?
        theta = lambda x: atan((r(x - 0.005) - r(x + 0.005)) / 0.01)

        # Находим координату, в которой угол наклона касательной будет больше 20 град, 
        # задаём точность расчёта и пересчитываем параметры относительного удлинения и относительного сужения кормовой части ??перепроверить это и возможно сделать аналитически
        if theta(self.L_cor) > radians(20):
            n = 0
            step_x = 0.01
            while theta(self.L - step_x * n) > radians(20):
                n += 1
                lambd_cor_fic = self.lambd * (self.L - step_x * n) / self.L
                eta_cor_fic = 2 * r(self.L - step_x * n) / self.D
            # Вызываем библиотечную функцию для определения c_xa0 кормовой части с параболической образующей,
            # имеющей lambd_cor и eta_cor, вычисленные по оставшейся части тела после отброса части, на которой
            # происходит отрыв потока (эта часть будет учтена в донном сопротивлении)
            return AeroBDSM.get_Cxcor_Par(M, lambd_cor, eta_cor)            
        else:
            # В данном случае по всей длине кормовой части не происходит отрыв потока
            # Вызываем библиотечную функцию для определения c_xa0 кормовой части с параболической образующей
            return AeroBDSM.get_Cxcor_Par(M, self.lambd, self.eta)        

    def coef_c_x_dn(self, M: float, is_active: bool) -> float:
        # Уравнение образующей кормовой части
        r = lambda x: self.D / 2 * (1 - self.eta) * (x / self.L)**2

        # Рассчитываем функцию угла наклона касательной к образующей  ??а нельзя это аналитически?
        theta = lambda x: atan((r(x - 0.005) - r(x + 0.005)) / 0.01)

        # Находим координату, в которой угол наклона касательной будет больше 20 град, 
        # задаём точность расчёта и пересчитываем параметры относительного удлинения и относительного сужения кормовой части ??перепроверить это и возможно сделать аналитически
        if theta(self.L_cor) > radians(20):
            n = 0
            step_x = 0.01
            while theta(self.L - step_x * n) > radians(20):
                n += 1
                lambd_cor_fic = self.lambd * (self.L - step_x * n) / self.L
                eta_cor_fic = 2 * r(self.L - step_x * n) / self.D

            # Тогда в качестве диаметра дна вместо фактического принимаем,
            # диаметр сечения, в котором произошёл отрыв потока
            D_dn_fic = eta_cor_fic * self.D
            S_dn_fic = pi * D_dn_fic**2 / 4

            # Вызываем библиотечную функцию для определения поправочного коэффициента, учитывающего сужение кормовой части
            # имеющей lambd_cor и eta_cor, вычисленные по оставшейся части тела после отброса части, на которой
            # происходит отрыв потока
            Sigma_eta = AeroBDSM.get_Sigma_eta(M, lambd_cor_fic, eta_cor_fic)
        
            # Вводим поправку на то, что площадь условного дна отличается от площади цилиндра
            if is_active == False:  #если сопло не работает учитываем всю площадь дна  
                return Sigma_eta * S_dn_fic / self.S_m
            else:                   #если сопло работает, то учитываем только площадь дна не занятую соплом
                return Sigma_eta * (S_dn_fic - self.S_a) / self.S_m
        else:
            # В данном случае по всей длине кормовой части не происходит отрыв потока
            # поэтому донное сопротивление рассчитывается по фактическому донному срезу кормовой части
            return super().coef_c_x_dn(M, is_active)
    
    def c_y_alpha(self, M: float) -> float:
        return super().c_y_alpha(M)

    def x_F_alpha(self, M: float) -> float:
        return super().x_F_alpha(M)
   
    def F_surf(self, x):
        super().F_surf(x)

        # Расчёт площади поверхности через интеграл по уравнению образующей линии
        return integrate.quad(lambda x: (self.D / 2 * (1 - self.eta) * (x / self.L)**2) * sqrt(
            1 + (self.D / 2 * (1 - self.eta) / self.L**2 * 2 * x)**2), 0, self.L)[0] * pi * 2
