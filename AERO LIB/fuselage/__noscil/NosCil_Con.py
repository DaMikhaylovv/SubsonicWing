from math import *

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM

from .NosCil import NosCil


class NosCil_Con(NosCil):
    '''Класс конической носовой части (конус + цилиндр)'''
    def __init__(self, *, lambd_nos: float, D: float, lambd_cil: float) -> None:
        # Расчёт основной геометрии (общей для любой носовой части)
        super().__init__(lambd_nos=lambd_nos, D=D, lambd_cil=lambd_cil)

        # Расчёт площади боковой проекции носовой части
        self.S_bok_nos = self.L_nos * D / 2

        # Расчёт объема носовой части
        self.W_nos = pi * D**2 * self.L_nos / 12

        # Расчёт полуугла раствора конуса
        self.theta_con = atan(1 / (2 * lambd_nos))

        # Расчёт координаты геометрического центра тяжести сечения носовой части
        self.x_cs_nos = 2 * D * lambd_nos / 3


    def c_xa0_p(self, M: float, H: float) -> float:       
        #Учитываем влияние пограничного слоя на конусе путём замены его фиктивным конусом с увеличенным углом раствора

        # Число Рейнольдса для конуса
        Re_con = Mach_to_v(M, H) * self.L_nos / atmo.nu(H)

        # Вызываем библиотечную функцию для определения числа Маха на поверхности конуса
        Ms_con = AeroBDSM.get_Ms_Con(M, self.theta_con)
        
        # Вычисляем дополнительный угол полураствора конуса 
        Delta_theta = 0.037 * (1 + 0.4 * Ms_con + 0.147 * Ms_con**2 - 0.006 * Ms_con**3) / Re_con**0.2

        # Угол полураствора фиктивного конуса
        theta_fic = self.theta_con + Delta_theta

        # Удлинение фиктивного конуса
        lambd_nos_fic = 1 /(2 * tan(theta_fic))

        # Вызываем библиотечную функцию для определения коэффициента сопротивления давления конической носовой части                
        return AeroBDSM.get_Cxnos_Con(M, lambd_nos_fic)  
    
    def Delta_c_x(self, M: float, alpha: float) -> float:
        '''
        Расчёт коэффициента дополнительной продольной силы, возникающей при ненулевых углах атаки,
        в результате перераспределения давления на носовой части

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                        
        Вывод:  Delta_c_x: float - коэффициент дополнительной продольной силы
        '''
        
        # Вызываем библиотечную функцию для определения коэффициента реализации дополнительной продольной силы на конической носовой части
        Sigma_c_xai_nos = AeroBDSM.get_Sigma_Cxinos(M, self.lambd_nos, 0)    

        return 2 * Sigma_c_xai_nos * sin(alpha)**2

    def c_y_alpha(self, M: float):
        # Вызываем библиотечную функцию для определения производной коэффициента нормальной силы по углу атаки
        return AeroBDSM.get_Cy_alpha_ConCil(M, self.lambd_nos, self.lambd_cil)

    def x_F_alpha(self, M: float) -> float:
        return super().x_F_alpha(M)  
        

    def lambd_nos_eqvCon(self, M: float) -> float:        
        # Носовая часть является конической, поэтому возвращаем удлинение носовой части
        return self.lambd_nos

    def kM(self, M: float) -> float:
        return super().kM(M)

    def F_surf(self, x: float) -> float:
        super().F_surf(x)

        if x == 0:
            return 0

        # Формула площади поверхности конуса в зависимости от координаты (начиная с носа ЛА)
        F_con = lambda x: pi * x**2 * self.D / (2 * self.L_nos * cos(self.theta_con))

        # Оцениваем положение введённой координаты и рассчитываем площадь боковой поверхности        
        if x <= self.L_nos:
            return F_con(x)
        return F_con(self.L_nos) + pi * self.D * (x - self.L_nos)