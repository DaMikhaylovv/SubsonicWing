from math import *
from numpy import sign

from libs.Atmosphere import atmo
from libs.handler import *
import libs.AeroBDSM as AeroBDSM

from .__noscil import *
from .__cor import *


class Fuselage:
    '''Класс фюзеляжа

    АДХ в связанной системе координат:
    c_x(M, alpha, H, is_active) - расчёт коэффициента продольной силы
    c_y(M, alpha) - расчёт коэффициента нормальной силы
    
    АДХ в скоростной системе координат:
    c_xa(M, alpha, H, is_active) - расчёт коэффициента силы лобового сопротивления
    c_ya(M, alpha, H, is_active) - расчёт коэффициента подъёмной силы

    x_d(M, alpha) - расчёт координаты центра давления

    Аэродинамические силы
    X(M, alpha, H, is_active) - продольная сила, Н
    Y(M, alpha, H) - нормальная сила, Н
    
    X_a(M, alpha, H, is_active) - сила лобового сопротивления, Н
    Y_a(M, alpha, H, is_active) - подъёмная сила, Н
    
    Вспомогательные АДХ    
    c_xa_dn(M, H, is_active) - расчёт коэффициента донного сопротивления
    c_xa0(M, H, is_active) - расчёт коэффициента силы лобового сопротивления при нулевом угле атаки
    c_xai(M, alpha) - расчёт коэффициента индуктивного сопротивления
    c_y_alpha(M) - расчёт производной по углу атаки коэффициента нормальной силы, 1/рад 
    x_F_alpha(M) - расчёт координаты фокуса по углу атаки

    Геометрия    
    F_surf(x) - расчёт площади поверхности без торцов   
    geometry() - вывод геометрических параметров
        
    '''
    def __init__(self, nosCil: NosCil, cor: Cor, h_s: float) -> None:
        # присваивание носовой и кормовой части
        self.nosCil = nosCil
        self.cor = cor

        # инициализация и расчёт геометрии фюзеляжа
        self.D = self.nosCil.D                                  # Диаметр миделя фюзеляжа
        self.L = self.nosCil.L + self.cor.L                     # Длина фюзеляжа
        self.S_m = self.nosCil.S_m                              # Площадь миделя
        self.lambd = self.L / self.D                            # Удлинение фюзеляжа

        self.h_s = h_s                                          # Высота неровностей (шероховатость) поверхности
        
        # вызываем библиотечную функцию для расчёта площади боковой проекции
        self.S_bok = AeroBDSM.get_S_bok(
            self.D,
            self.nosCil.S_bok_nos,
            self.nosCil.lambd_cil,
            self.cor.lambd,
            self.cor.eta,
            0,                  #несущих поверхностей и областей затенения нет
            0,
            self.cor.L).Value
              
        self.F_s = self.F_surf(self.L)           # Передача площади поверхности фюзеляжа (без торцов)
    
    # @checkSavedValue
    def c_y_alpha(self, M: float) -> float:
        '''
        Расчёт производной по углу атаки коэффициента нормальной силы фюзеляжа, 1/рад
        
        Ввод:   M: float - число Маха
                        
        Вывод:  c_y_alpha: float - производная по углу атаки коэффициента нормальной силы фюзеляжа, 1/рад
        '''
        return self.nosCil.c_y_alpha(M) + self.cor.c_y_alpha(M)
    
    # @checkSavedValue
    def __kappa_a(self, M: float, alpha: float) -> float:
        '''
        Эмпирическая поправка, учитывающая замедление роста "безотрывной" составляющей 
        нормальной силы фюзеляжа по мере увеличения угла атаки (альфа в градусах а можно просто 2pi?)
        '''
        #return 1 - 0.45 * ( 1 - exp( -0.06 * M**2 ) * ( 1 - exp( -2.1885 * pi * abs(alpha) ) ) )
        return 1 - 0.45 * ( 1 - exp( -0.06 * M**2 ) * ( 1 - exp( -2 * pi * abs(alpha) ) ) )
        #return 1 - 0.45 * ( 1 - exp( -0.06 * M**2 ) * ( 1 - exp( -0.12 * abs(alpha) ) ) )

    # @checkSavedValue
    def c_y(self, M: float, alpha: float) -> float:
        '''
        Расчёт коэффициента нормальной силы фюзеляжа

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад

        Вывод:  c_y - коэффициент нормальной силы
        '''
        
        #?
        # S_bok будем считать при b_b, l_hv_1 = 0
        # Считаем в интерференции затеняемую площадь для внесения поправок
        # S_bok_z = S_bok(b_b, l_hv_1 = 0) - S_bok(b_b, l_hv_1)
        # Для c_y1_f следует вычесть err = 4 * S_bok_z * c_x_cil * (sin(alpha))**2 * sign(alpha) / (pi * D**2)
        # Для c_y_f вычитаем err * cos(alpha)        
        # L_hv_1 = g['L_hv'] - g['b_b'] / 2 - на данном этапе зануляем

        #Вызываем библиотечную функцию для определения коэффициента сопротивления цилиндра при обтекании по нормали к его оси
        c_x_cil_N = AeroBDSM.get_Cx_Cil_N(M, abs(alpha))
        
        return self.c_y_alpha(M) * self.__kappa_a(M, alpha) * sin(alpha) * cos(alpha) + \
            c_x_cil_N * sin(alpha)**2 * sign(alpha) * self.S_bok / self.S_m

    # @checkSavedValue
    def c_xa0_p(self, M: float, H: float) -> float:
        '''
        Расчёт коэффициента лобового сопротивления давления при нулевом угле атаки

        Ввод:   M: float - число Маха
                H: float - высота полета, м    

        Вывод:  c_xa0_p - коэффициент лобового сопротивления давления при нулевом угле атаки
        '''       
        
        return self.nosCil.c_xa0_p(M, H) + self.cor.c_xa0_p(M)
    
    # @checkSavedValue
    def x_t(self, M: float, H: float, T_s: float = None) -> float:
        '''
        Расчёт координаты точки перехода ламинарного пограничного слоя в турбулентный
        
        Ввод:   M: float - число Маха
                H: float - высота полета, м
                T_s: float - средняя температура поверхности фюзеляжа, К

        Вывод:  x_t: float - координата точки перехода ламинарного пограничного слоя в турбулентный (от носа), м
        '''         
       
        # Относительная высота бугорков шероховатости поверхности
        hh = self.h_s / self.L

        # Число Рейнольдса для фюзеляжа
        Re_f = Mach_to_v(M, H) * self.L / atmo.nu(H)

        # Вызываем библиотечную функцию для определения критического числа Рейнольдса
        # при отсутствии теплопередачи и нулевом градиенте давления
        Re_t_0_1 = AeroBDSM.get_Re_t0_1(M, Re_f, hh)

        #Это не понятный кусок, надо ли это? Получается довольно странно, нужны исследования
        # Если носовая часть имеет выпуклую образующую, то при М>1 Re_t_0 следует
        # увеличить на 30..50%. Для этого введём поправочный множитель kW, зависящий от
        # объёма носовой части. Для объёма конической носовой части kW=1 (наименьший объём,
        # прямая образующая), для сферической kW=1.5  Это надо сделать
        #kW = 1.3
        #Будем вводить эту поправку линейно по М, от 0.8 (может Мкр фюзеляжа),
        #до М = 1  
        #if M <= 0.8:
        #    sigma_W = 1.0
        #elif 0.8 < M < 1.0:
        #    sigma_W = 1.0 * ((1.0 - M) / (1.0 - 0.8)) + kW * (1 - (1.0 - M) / (1.0 - 0.8))
        #else:
        #    sigma_W = kW
        

        sigma_W = 1
        
        Re_t_0 = Re_t_0_1 * sigma_W
       
        if T_s == None:
            # Средняя температура носовой части считается равной температуре восстановления (теплопередача отсутствует)
            # В этом случае, относительная температура стенки равна 1 (отношение температуры стенки к температуре восстановления)
            TT_s = 1
        else:
            # Задана известная средняя температура поверхности в носовой части фюзеляжа
            # Рассчитываем температуру восстановления            
            # r = 0.845..0.88 - коэффициент восстановления температуры для ламинарного и турбулентного слоя соответственно
            # возьмём усреднённое значение r = 0.865  а какой надо брать то?
            T_r = atmo.T(H) * (1 + (atmo.k - 1) / 2 * 0.865 * M**2)
            # Относительная температура стенки
            TT_s = self.T_s / T_r

        # Вызываем библиотечную функцию для определения поправочного множителя,
        # учитывающего влияние температуры поверхности тела на критическое число Рейнольдса
        Sigma_T = AeroBDSM.get_Sigma_T(M, TT_s).Value
        
        # Критическое число Рейнольдса
        Re_t = Re_t_0 * Sigma_T

        # Вычислим координату точки перехода ламинарного пограничного слоя в турбулентный
        x_t = Re_t * atmo.nu(H) / Mach_to_v(M, H)
         
        return x_t    
    
    # @checkSavedValue
    def F_surf(self, x):
        '''
        Расчёт площади поверхности (без торцов)
        
        Ввод:   x: float - координата границы области расчёта площади поверхности, м
        
        Вывод:  F_surf: float - площадь поверхности, м^2
        '''
        assert all([x >= 0, x <= self.L]), \
            RuntimeError(f'Используемая координата выходит за пределы фюзеляжа, значение координаты: {x}')
        if x <= self.nosCil.L:
            return self.nosCil.F_surf(x)
        return self.nosCil.F_surf(self.nosCil.L) + self.cor.F_surf(round(x - self.nosCil.L, 8)) #? зачем тут self.nosCil.L, 8
    
    # @checkSavedValue
    def c_xa0_f(self, M: float, H: float, x_t: float = None, T_s: float = None) -> float:
        '''
        Расчёт коэффициента лобового сопротивления трения при нулевом угле атаки

        Ввод:   M: float - число Маха
                H: float - высота полета, м
                x_t: float - координата точки перехода ламинарного пограничного слоя в турбулентный (от носа), м
                T_s: float - средняя температура поверхности фюзеляжа, К

        Вывод:  c_xa0_f - коэффициент лобового сопротивления трения при нулевом угле атаки
        '''       
        #??Что насчёт поправки к тернию на конус?
        # Число Рейнольдса фюзеляжа
        Re_f = Mach_to_v(M, H) * self.L / atmo.nu(H)

        if x_t == None:
            #Если координата точки перехода ЛПС в ТПС не задана, то определяем её по критическому числу Рейнольдса
            x_t = self.x_t(M, H, T_s)

        # Относительная координата точки перехода ламинарного пограничного слоя в турбулентный
        xx_t = self.F_surf(x_t) / self.F_s

        if T_s == None:
            # Средняя температура носовой части считается равной температуре восстановления (теплопередача отсутствует)
            # Вызываем библиотечную функцию для определения коэффициента трения плоской пластинки при М = 0
            cf_M0 = AeroBDSM.get_Cf_M0(Re_f, xx_t)
            
            # Вызываем библиотечную функцию для определения поправочного множителя, учитывающего влияние числа Маха на коэффициент трения плоской пластинки
            sigma_M = AeroBDSM.get_Sigma_M(M, xx_t)
            
            # Коэффициент трения плоской пластинки с учётом влияния числа Маха
            cf = cf_M0 * sigma_M
        else:
            # Задана известная средняя температура поверхности фюзеляжа
            # Рассчитываем температуру восстановления            
            # r = 0.845..0.88 - коэффициент восстановления температуры для ламинарного и турбулентного слоя соответственно
            # возьмём усреднённое значение r = 0.865  ?? а какой надо брать то?
            T_r = atmo.T(H) * (1 + (atmo.k - 1) / 2 * 0.865 * M**2)
            # Рассчитываем определяющую температуру 
            T_zv = (atmo.T(H) + T_s) / 2 + 0.22 * (T_r - T_s)
            # Число Рейнольдса, соответствующее определяющей температуре 
            Re_zv = Re_f / (T_zv / atmo.T(H)**1.76)
            # Вызываем библиотечную функцию для определения коэффициента трения плоской пластинки при М = 0
            cf_M0 = AeroBDSM.get_Cf_M0(Re_zv, xx_t)
            # Вызываем библиотечную функцию для определения поправочного множителя, учитывающего влияние числа Маха на коэффициент трения плоской пластинки
            sigma_M = AeroBDSM.get_Sigma_M(M, xx_t) #?? В ЛиЧ не написано что здесь надо поправку на М
            
            # Коэффициент трения плоской пластинки с учётом влияния теплопередачи и числа Маха
            cf = cf_M0 * atmo.T(H) / T_zv * sigma_M
                   
        # Вычисляем коэффициент сопротивления трения фюзеляжа  
        return cf * self.F_s / self.S_m

    # @checkSavedValue
    def c_xa_dn(self, M: float, H: float, is_active: bool, x_t: float = None, T_s: float = None) -> float:
        '''
        Расчёт коэффициента донного сопротивления

        Ввод:   M: float - число Маха
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)
                x_t: float - координата точки перехода ламинарного пограничного слоя в турбулентный (от носа), м
                T_s: float - средняя температура поверхности фюзеляжа, К

        Вывод:  c_xa_dn - коэффициент донного сопротивления
        '''

        # Расчет параметров атмосферы для заданной высоты
        nu = atmo.nu(H)
        v = Mach_to_v(M, H)
        a = atmo.a(H)

        # Расчёт числа Рейнольдса при числе Маха М = 0.6
        Re_06M = 0.6 * a * self.L / nu

        # Расчёт числа Рейнольдса
        Re_f = v * self.L / nu

        if x_t == None:
            #Если координата точки перехода ЛПС в ТПС не задана, то определяем её по критическому числу Рейнольдса
            x_t = self.x_t(M, H, T_s)

        # Относительная координата точки перехода ламинарного пограничного слоя в турбулентный
        xx_t = self.F_surf(x_t) / self.F_s
        
        # c_x_dn будет считаться при c_prof = 0 ?
        # т.к. с отличным от 0 c_prof значение c_x_dn будет выше, следует в разделе интерференции сделать 
        # прибавку к зависимым результатам класса FUS
        # Посчитаем delta_c_x_dn = c_x_dn(c_prof /= 0) - c_x_dn(c_prof = 0)
        # c_x0_f станет алгебраически больше на delta_c_x_dn, следовательно 
        # из c_y_f = c_y1_f * cos(alpha) - c_x0_f * sin(alpha) следует вычесть delta_c_x_dn * sin(alpha)
        
        # Вызываем библиотечную функцию для определения коэффициента донного сопротивления цилиндра
        c_xa_dn_cil = AeroBDSM.get_Cxdon_Cil(M, 0, Re_06M, Re_f, xx_t, self.L / self.D).Value

        # Вычисляем коэффициент донного сопротивления с учётом сопла и сужения кормовой части
        return c_xa_dn_cil * self.cor.coef_c_x_dn(M, is_active)

    # @checkSavedValue
    def c_xa0(self, M: float, H: float, is_active: bool, x_t: float = None, T_s: float = None) -> float:
        '''
        Расчёт коэффициента лобового сопротивления давления при нулевом угле атаки

        Ввод:   M: float - число Маха
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  c_xa0 - коэффициент лобового сопротивления при нулевом угле атаки
        '''       
        
        return self.c_xa0_p(M, H) + self.c_xa0_f(M, H, x_t, T_s) + self.c_xa_dn(M, H, is_active, x_t, T_s)
    
    # @checkSavedValue
    def c_xai(self, M: float, alpha: float) -> float:
        '''
        Расчёт коэффициента индуктивного сопротивления

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад

        Вывод:  c_xi - коэффициент индуктивного сопротивления
        '''

        #Вычисляем коэффициент дополнительной продольной силы носовой части,
        #возникающей при углах атаки
        Delta_c_x = self.nosCil.Delta_c_x(M, alpha)

        return self.c_y(M, alpha) * sin(alpha) + Delta_c_x * cos(alpha)

    # @checkSavedValue
    def c_xa(self, M: float, alpha: float, H: float, is_active: bool, x_t: float = None, T_s: float = None) -> float:
        '''
        Расчёт коэффициента силы лобового сопротивления

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  c_xa - коэффициент силы лобового сопротивления
        '''
        # При расчёте интерференции стоит прибавить c_xi_1 * S_r_kon / S_m * k_tau1 
        # и 1.05 * c_x0_1 * k_tau1 * S_r_kon / S_m ?????
        
        return self.c_xa0(M, H, is_active, x_t, T_s) + self.c_xai(M, alpha)

    # @checkSavedValue
    def c_x(self, M: float, alpha: float, H: float, is_active: bool, x_t: float = None, T_s: float = None) -> float:
        '''
        Расчёт коэффициента продольной силы (пересчёт в связанную СК)

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  c_x - коэффициент продольной силы
        '''
        return (self.c_xa(M, alpha, H, is_active, x_t, T_s) - self.c_y(M, alpha) * sin(alpha)) / cos(alpha)

    # @checkSavedValue
    def c_ya(self, M: float, alpha: float, H: float, is_active: bool, x_t: float = None, T_s: float = None) -> float:
        '''
        Расчёт коэффициента подъёмной силы (пересчёт в скоростную СК)

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  c_ya - коэффициент подъёмной силы
        '''
        return (self.c_y(M, alpha) - self.c_xa(M, alpha, H, is_active, x_t, T_s) * sin(alpha) ) / cos(alpha)        

    # @checkSavedValue
    def x_F_alpha(self, M: float) -> float:
        '''
        Расчёт координаты фокуса по углу атаки

        Ввод:   M: float - число Маха
                        
        Вывод:  x_F_alpha: float - координата фокуса по углу атаки (от носа), м
        '''
        return (self.nosCil.x_F_alpha(M) * self.nosCil.c_y_alpha(M) +
                (self.cor.x_F_alpha(M) + self.nosCil.L_nos + self.nosCil.L_cil) * self.cor.c_y_alpha(M)) / self.c_y_alpha(M)

    # @checkSavedValue
    def x_d(self, M: float, alpha: float) -> float:
        '''
        Расчёт координаты центра давления

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                        
        Вывод:  x_d: float - координата центра давления (от носа), м
        '''
        
        # x_c_S_bok здесь считается при b_b = 0 и L_hv = 0, ?
        # будем считать погрешность delta_x_c_S_bok и прибавлять в разделе интерференции,
        # домножая на 4 * self.S_bok(M) * self.c_x_cil * (sin(alpha))**2 * sign(alpha) / (pi * self.D**2)
        # учитывая поправку на S_bok
        
        # При нулевой нормальной силе центр давления отсутствует, условно принимаем что центр давления совпадает с фокусом
        if self.c_y(M, alpha) == 0:
            return self.x_F_alpha(M)

        # Вызываем библиотечную функцию для расчёта центра тяжести площади боковой проекции фюзеляжа
        x_c_S_bok = AeroBDSM.get_x_c_pl_F(
            self.D,
            self.nosCil.S_bok_nos,
            self.nosCil.x_cs_nos,
            self.nosCil.lambd_nos,
            self.nosCil.lambd_cil,
            self.cor.lambd,
            self.cor.eta,
            M,
            0,
            self.cor.L).Value
        
        #Вызываем библиотечную функцию для определения коэффициента сопротивления цилиндра при обтекании по нормали к его оси
        c_x_cil_N = AeroBDSM.get_Cx_Cil_N(M, abs(alpha))

        # Предполагая, что линейная по углу атаки составляющая нормальной силы приложена в фокусе,
        # а нелинейная - в центре тяжести площади боковой проекции, рассчитываем центр давления
        # суммарной силы таким образом, чтобы создаваемый ею момент был равен сумме моментов от двух составляющих        
        return (self.c_y_alpha(M) * self.__kappa_a(M, alpha) * sin(alpha) * cos(alpha) * self.x_F_alpha(M) +
                 c_x_cil_N * sin(alpha)**2 * sign(alpha) * x_c_S_bok * self.S_bok / self.S_m ) / self.c_y(M, alpha)

    @checkSavedValue
    def m_z_omega_z(self, M: float, xx_0: float) -> float:
        '''
        Расчет коэффициента продольного демпфирующего момента

        Ввод:   M: float - число Маха
                xx_0: float - относительная координата центра вращения в долях длины фюзеляжа (от носа) 
        Вывод:  m_z_omega_z: float - коэффициент продольного демпфирующего момента, с/рад
        '''

        # Координата центра тяжести объёма фюзеляжа (? эта формула для конус+цилиндр, надо сделать более точно для каждого типа носа и кормы)
        x_cV = self.L * (2 * self.lambd**2 - self.nosCil.lambd_nos**2) / (4 * self.lambd * (self.lambd - 2 / 3 * self.nosCil.lambd_nos))

        m_z_omega_z = - 2 * (1 - xx_0 + xx_0**2 - x_cV / self.L)

        # Результат
        return m_z_omega_z
    
    @checkSavedValue
    def m_z(self, M: float, alpha: float, xx_0: float, omega_z: float = 0, H: float = 0) -> float:
        '''
        Расчет коэффициента продольного момента

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                xx_0: float - относительная координата точки, относительно которой определяется момент (в долях длины фюзеляжа, от носа) 
                omega_z: float - угловая скорость вращения фюзеляжа вокруг оси Оz, рад/с
                H: float - высота полета, м
        Вывод:  m_z: float - коэффициент продольного момента
        '''
        
        # Координата центра давления
        x_d = self.x_d(M, alpha)

        # Координата точки, относительно которой определяется момент
        x_0 = self.L * xx_0

        # Коэффициент продольного демпфирующего момента
        m_z_omega_z = self.m_z_omega_z(M, xx_0)

        # Безразмерная угловая скорость
        omega_omega_z = omega_z * self.L / Mach_to_v(M, H)

        # Коэффициент продольного момента
        m_z = self.c_y(M, alpha) * (x_0 - x_d) / self.L + m_z_omega_z * omega_omega_z

        # Результат
        return m_z    
   
    def geometry(self):
        '''
        Функция вывода геометрии. ??надо чтобы эта ф-ция вызывала ф-ции подклассов и печатали все геом. параметры

        Вывод: Dict[str, float], где ключи -- названия геометрического параметра
        '''
        attrlist = [
            'D', 'lambd_nos', 'lambd_cil', 'r_sph', 'S_bok_nos', 'W_nos', 'x_cs_nos', 'L_cil', 'L_nos', 'S_m', 'L', 'rr_sph',
            'L_cor', 'lambd_cor', 'W_cor', 'eta_cor', 'D_dn', 'D_a']
        if hasattr(self, 'theta_nos_con'):
            attrlist.append('theta_nos_con')
        return {attrname:self.__getattribute__(attrname) for attrname in attrlist}
      
    # @checkSavedValue
    def X(self, M: float, alpha: float, H: float, is_active: bool) -> float:
        '''
        Расчёт продольной силы

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  X: float - продольная сила, Н
        '''
        return self.c_x(M, alpha, H, is_active) * self.S_m * atmo.rho(H) * Mach_to_v(M, H)**2 / 2
  
    # @checkSavedValue
    def Y(self, M: float, alpha: float, H: float) -> float:
        '''
        Расчёт подъёмной силы в связанной системе координат

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м

        Вывод:  Y: float - нормальная сила, Н
        '''
        return self.c_y(M, alpha) * self.S_m * atmo.rho(H) * Mach_to_v(M, H)**2 / 2
    
    # @checkSavedValue
    def X_a(self, M: float, alpha: float, H: float, is_active: bool) -> float:
        '''
        Расчёт силы лобового сопротивления

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  X_a: float - сила лобового сопротивления, Н
        '''
        return self.c_xa(M, alpha, H, is_active) * self.S_m * atmo.rho(H) * Mach_to_v(M, H)**2 / 2

    # @checkSavedValue
    def Y_a(self, M: float, alpha: float, H: float, is_active: bool) -> float:
        '''
        Расчёт подъёмной силы

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                H: float - высота полета, м
                is_active: bool - режим работы двигателя: активный (True), пассивный (False)

        Вывод:  Y_a: float - подъёмная сила, Н
        '''
        return self.c_ya(M, alpha, H, is_active) * self.S_m * atmo.rho(H) * Mach_to_v(M, H)**2  / 2
    
    @checkSavedValue
    def M_z(self, M: float, alpha: float, xx_0: float, omega_z: float = 0, H: float = 0) -> float:
        '''
        Расчет продольного момента

        Ввод:   M: float - число Маха
                alpha: float - угол атаки, рад
                xx_0: float - относительная координата точки, относительно которой определяется момент (в долях длины фюзеляжа, от носа) 
                omega_z: float - угловая скорость вращения, рад/с
                H: float - высота полета, м
        Вывод:  M_z: float - продольный момент, Н*м
        '''

        # Коэффициент продольного момента
        m_z = self.m_z(M, alpha, xx_0, omega_z, H)

        # Продольный момент
        M_z = m_z * self.S * self.L * q(M, H)

        # Результат
        return M_z