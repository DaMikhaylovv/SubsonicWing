from math import *

import scipy.integrate as integrate
from scipy.optimize import bisect
import numpy as np

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM

from .NosCil import NosCil


class NosCil_SphOgi(NosCil):
    '''Класс оживальной носовой части с сферическим затуплением (сфера + оживало + цилиндр)'''
    def __init__(self, *, lambd_nos: float, D: float, lambd_cil: float, r_sph: float) -> None:
        # Расчёт основной геометрии (общей для любой носовой части)
        super().__init__(lambd_nos=lambd_nos, D=D, lambd_cil=lambd_cil, r_sph=r_sph)

        # Вспомогательная переменная
        a = self.L_nos - r_sph

        # Угол сектора окружности-оживала, занимаемого оживальной частью носа (вспомогательный угол)
        angle = np.arccos( (a**2 - (D / 2 - r_sph)**2) / (a**2 + (D / 2 - r_sph)**2) )

        # Координата перехода сфера-оживало
        x1 = r_sph * (1 - sin(angle))

        # Радиус оживала
        R = r_sph + a / sin(angle)

        ### Нахождение площади продольного сечения носовой части

        # Площадь сферической части сечения
        S_bok_sph = integrate.quad(lambda x: sqrt(r_sph**2 - (x - r_sph)**2), 0, x1)[0] * 2

        # Площадь оживальной части сечения
        S_bok_ogi = integrate.quad(lambda x: sqrt(R**2 - (x - self.L_nos)**2) - R + D / 2, x1, self.L_nos)[0] * 2

        ### Нахождение объема носовой части

        # Объем сферической части
        W_sph = integrate.quad(lambda x: r_sph**2 - (x - r_sph)**2, 0, x1)[0] * pi

        # Объем оживальной части
        W_ogi = integrate.quad(lambda x: (sqrt(R**2 - (x - self.L_nos)**2) - R + D / 2)**2, x1, self.L_nos)[0] * pi

        # Сохранение результата расчета
        self.S_bok_nos = S_bok_sph + S_bok_ogi
        self.W_nos = W_sph + W_ogi

        # Функция зависимости радиуса носовой части от координаты Х
        def r(x):
            if x <= x1:
                return sqrt(r_sph**2 - (x - r_sph)**2)
            return D / 2 - R + sqrt( R**2 - (x - self.L_nos)**2 )

        # Расчёт координаты геометрического центра тяжести сечения носовой части
        self._x_cs_nos(r)

    # @checkSavedValue
    def __c_xa0_p_sphOgi(self, M: float) -> float:
        # Расчёт удлинения незатупленного оживала
        lambd_nos_ogi = (self.lambd_nos - self.rr_sph / 2) / sqrt(1 - self.rr_sph)

        # Вызываем библиотечную функцию для определения коэффициента сопротивления давления параболической носовой части
        c_xa_ogi = AeroBDSM.get_Cxnos_Par(M, lambd_nos_ogi)

        # Вызываем библиотечную функцию для определения коэффициента сопротивления давления сферической носовой части
        c_xa_sph = AeroBDSM.get_Cxnos_Ell(M, 0.5)

        # Угол наклона образующей в точке сопряжения сферы с оживалом
        theta = atan(sqrt(1 - self.rr_sph) / lambd_nos_ogi)

        # Расчёт коэффициента лобового сопротивления носовой части сфера+оживало
        return c_xa_ogi * (1 - self.rr_sph**2 * cos(theta)**2 * (3.1 - 1.4 * self.rr_sph * cos(theta) - 0.7 * self.rr_sph**2 * cos(theta)**2)) + c_xa_sph * self.rr_sph**2


    # @checkSavedValue
    def c_xa0_p(self, M: float, H: float) -> float:      
        # Расчёт коэффициента сопротивления давления носовой части сфера+оживало (без учёта пограничного слоя)
        c_xa_sphOgi = self.__c_xa0_p_sphOgi(M)

        ##Учитываем влияние пограничного слоя путём вычисления поправки на основе эквивалентной конической носовой части (или не надо? или надо но по другому??)

        # Определяем удлинение эквивалентной конической носовой части
        lambd_nos_eq_con = self.lambd_nos_eqvCon(M)
       
        # Вычисляем длину эквивалентной конической носовой части
        L_eq_con = lambd_nos_eq_con * self.D
        
        # Вычисляем угол полураствора эквивалентного конуса
        theta_eq_con = atan(self.D / (2 * L_eq_con))

        # Число Рейнольдса для конуса
        Re_con = Mach_to_v(M, H) * L_eq_con / atmo.nu(H)

        # Вызываем библиотечную функцию для определения числа Маха на поверхности конуса
        Ms_con = AeroBDSM.get_Ms_Con(M, theta_eq_con)
        
        # Вычисляем дополнительный угол полураствора конуса 
        Delta_theta = 0.037 * (1 + 0.4 * Ms_con + 0.147 * Ms_con**2 - 0.006 * Ms_con**3) / Re_con**0.2

        # Угол полураствора фиктивного конуса
        theta_fic = theta_eq_con + Delta_theta

        # Удлинение фиктивного конуса
        lambd_nos_fic = 1 /(2 * tan(theta_fic))
             
       
        # Вызываем библиотечную функцию для определения коэффициента сопротивления давления эквивалентной конической носовой части ??оно ведь совпадает с c_xa_sphOgi?? нет, потому что могли упереться в самый тупой конус
        c_xa_nos_eq_con = AeroBDSM.get_Cxnos_Con(M, lambd_nos_eq_con)

        # Вызываем библиотечную функцию для определения коэффициента сопротивления давления фиктивной конической носовой части
        c_xa_nos_fic_con = AeroBDSM.get_Cxnos_Con(M, lambd_nos_fic) 

        # Вычисляем коэффициент, учитывающий влияние пограничного слоя
        k_delta = c_xa_nos_fic_con / c_xa_nos_eq_con

        #print(f'''{M = }\n{k_delta = }\n{c_xa_nos_fic_con = }\n{c_xa_nos_eq_con = }\n{c_xa_sphOgi = }\n{lambd_nos_fic = }\n{lambd_nos_eq_con = }''')

        # Вычисляем коэффициент сопротивления давления носовой части сфера+оживало с учётом влияния пограничного слоя
        return c_xa_sphOgi * k_delta   
   
    def Delta_c_x(self, M: float, alpha: float) -> float:
        '''
        Расчёт коэффициента дополнительной продольной силы, возникающей при ненулевых углах атаки,
        в результате перераспределения давления на носовой части

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                        
        Вывод:  Delta_c_x: float - коэффициент дополнительной продольной силы
        '''
        
        # Вызываем библиотечную функцию для определения коэффициента реализации дополнительной продольной силы на оживальной носовой части
        Sigma_c_xai_nos = AeroBDSM.get_Sigma_Cxinos(M, self.lambd_nos, 1)

        return 2 * Sigma_c_xai_nos * sin(alpha)**2

    # @checkSavedValue
    def c_y_alpha(self, M: float) -> float:
        # Расчёт удлинения незатупленного оживала
        lambd_nos_ogi = (self.lambd_nos - self.rr_sph / 2) / sqrt(1 - self.rr_sph)

        # Вызываем библиотечную функцию для определения коэффициента нормальной силы оживало+цилиндр
        c_y_alpha_ogi = AeroBDSM.get_Cy_alpha_OgiCil(M, lambd_nos_ogi, self.lambd_cil)

        # Вызываем библиотечную функцию для определения коэффициента нормальной силы сфера+цилиндр
        c_y_alpha_sph = AeroBDSM.get_Cy_alpha_SphCil(M, 0.5, self.lambd_cil)
        
        return c_y_alpha_ogi * (1 - self.rr_sph**2) + c_y_alpha_sph * self.rr_sph**2

    def x_F_alpha(self, M: float) -> float:
        return super().x_F_alpha(M)

    # @checkSavedValue
    def lambd_nos_eqvCon(self, M: float) -> float:   
        # Расчёт коэффициента сопротивления давления носовой части сфера+оживало (без учёта пограничного слоя)
        c_xa_sphOgi = self.__c_xa0_p_sphOgi(M)

        # Расчёт коэффициента лобового сопротивления эквивалентной носовой части
        c_x_nos_eq_con = AeroBDSM.get_Cxnos_Con(M, 1.5)

        borders = (c_x_nos_eq_con.InputComplexes[1].Min,
                   c_x_nos_eq_con.InputComplexes[1].Max)
        
        func = lambda x: AeroBDSM.get_Cxnos_Con(M, x).Value - c_xa_sphOgi

        if func(borders[0]) < 0:
            lambd_nos_eq_con = borders[0]
        elif func(borders[1]) > 0:
            lambd_nos_eq_con = borders[1]
        else:
            lambd_nos_eq_con = bisect(func, *borders, disp=True, xtol=1e-4)

        # Возвращение результата        
        return lambd_nos_eq_con

    def kM(self, M: float) -> float:
        return super().kM(M)

    # @checkSavedValue
    def F_surf(self, x: float) -> float:
        super().F_surf(x)

        if x == 0:
            return 0

        # Вспомогательные переменные
        r_sph = self.r_sph
        D = self.D
        a = self.L_nos - r_sph
        L_nos = self.L_nos

        # Угол сектора окружности-оживала, занимаемого оживальной частью носа (вспомогательный угол)
        tmp = (D / 2 - r_sph)
        angle = np.arccos( (a * a - tmp * tmp) / (a * a + tmp * tmp) )

        # Координата перехода сфера-оживало (вспомогательная переменная)
        x1 = r_sph * (1 - cos(angle))

        # Радиус оживала
        R = r_sph + a / sin(angle)

        # Формула расчета боковой площади сферического затупления
        F_sph = lambda x: 2 * pi * r_sph * x
        
        # Площадь оживальной части сечения       
        F_ogi = lambda x: integrate.quad(lambda x: R * (1 - (R - D / 2) / sqrt((R - x + L_nos) * (R + x - L_nos))), x1, x)[0] * 2 * pi

        # Оцениваем положение введённой координаты и расчитываем площадь боковой поверхности
        if x <= x1:
            return F_sph(x)
        if x <= self.L_nos:
            return F_sph(x1) + F_ogi(x)
        return F_sph(x1) + F_ogi(self.L_nos) + pi * self.D * (x - self.L_nos)
