from math import *
import scipy.integrate as integrate
from scipy.optimize import bisect

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM

from .NosCil import NosCil


class NosCil_Ell(NosCil):
    '''Класс эллиптической носовой части (эллипсоид + цилиндр)'''
    def __init__(self, *, lambd_nos: float, D: float, lambd_cil: float) -> None:
        # Расчёт основной геометрии (общей для любой носовой части)
        super().__init__(lambd_nos=lambd_nos, D=D, lambd_cil=lambd_cil)

        # Расчёт площади продольного сечения носовой части
        self.S_bok_nos = integrate.quad(
                lambda x: sqrt((D / 2)**2 - (x - self.L_nos)**2 * (D / 2 / self.L_nos)**2), 0, self.L_nos)[0] * 2

        # Расчёт объема носовой части
        self.W_nos = integrate.quad(
                lambda x: sqrt((D / 2)**2 - (x - self.L_nos)**2 * (D / 2 / self.L_nos)**2)**2, 0, self.L_nos)[0] * pi

        # Расчёт координаты геометрического центра тяжести сечения носовой части
        self.x_cs_nos = (3 * pi - 4) * D * lambd_nos / (3 * pi)

    # @checkSavedValue
    def c_xa0_p(self, M: float, H: float) -> float:
        # Вызываем библиотечную функцию для определения коэффициента сопротивления давления эллиптической носовой части
        c_xa_Ell = AeroBDSM.get_Cxnos_Ell(M, self.lambd_nos)

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

        #print(f'''{M = }\n{k_delta = }\n{c_xa_nos_fic_con = }\n{c_xa_nos_eq_con = }\n{c_xa_Ell = }\n{lambd_nos_fic = }\n{lambd_nos_eq_con = }''')

        # Вычисляем коэффициент сопротивления давления эллиптической носовой части с учётом влияния пограничного слоя
        return c_xa_Ell * k_delta   

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
        # Вызываем библиотечную функцию для определения коэффициента нормальной силы эллиптической носовой части
        return AeroBDSM.get_Cy_alpha_SphCil(M, self.lambd_nos, self.lambd_cil)  

    def x_F_alpha(self, M: float) -> float:
        return super().x_F_alpha(M)
      
    def lambd_nos_eqvCon(self, M: float) -> float:
        
        # Расчёт коэффициента лобового сопротивления носовой части (без учёта дальнейших поправок)
        c_x_nos = AeroBDSM.get_Cxnos_Ell(M, self.lambd_nos)

        # Расчёт коэффициента лобового сопротивления эквивалентной носовой части (немножко чуда)
        c_x_nos_eq_con = AeroBDSM.get_Cxnos_Con(M, 1.5)

        borders = (c_x_nos_eq_con.InputComplexes[1].Min,
                   c_x_nos_eq_con.InputComplexes[1].Max)
        
        func = lambda x: AeroBDSM.get_Cxnos_Con(M, x).Value - c_x_nos

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

    def F_surf(self, x: float) -> float:
        super().F_surf(x)

        if x == 0:
            return 0

        # Вспомогательные переменные
        a = self.L_nos
        b = self.D / 2

        # Формула площади поверхности эллиптической носовой части
        F_ell = lambda x: integrate.quad(
            lambda x: sqrt(b**2 - (x - a)**2 * (b / a)**2) * sqrt(1 + (b * (a - x) / (a**2 * sqrt(1 - (x - a)**2 / a**2)))**2), 0, x)[0] * 2 * pi

        # Оцениваем положение введённой координаты и рассчитываем площадь боковой поверхности        
        if x <= self.L_nos:
            return F_ell(x)        
        return F_ell(self.L_nos) + pi * self.D * (x - self.L_nos)