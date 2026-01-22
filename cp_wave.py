import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class ThreeWaveMixing:
    def __init__(self, 
                 lambda1_nm, lambda2_nm, 
                 n1, n2, n3, 
                 deff_pm_V, 
                 L_mm,
                 process_type='SFG'):
        """
        初始化非线性晶体参数
        :param lambda1_nm: 信号光波长 (nm)
        :param lambda2_nm: 闲频光波长 (nm)
        :param n1, n2, n3: 三个波长对应的折射率
        :param deff_pm_V: 有效非线性系数 (pm/V)
        :param L_mm: 晶体长度 (mm)
        :param process_type: 非线性过程类型 ('SFG' 或 'SHG')
        """
        # --- 1. 物理常量 (SI Units) ---
        self.c = 2.99792458e8       # 光速 (m/s)
        self.epsilon0 = 8.854e-12   # 真空介电常数 (F/m)
        
        # --- 2. 参数转换 (Engineering -> SI) ---
        self.L = L_mm * 1e-3        # mm -> m
        self.deff = deff_pm_V * 1e-12 # pm/V -> m/V
        self.process_type = process_type
        
        # 波长与频率
        self.lam1 = lambda1_nm * 1e-9
        self.lam2 = lambda2_nm * 1e-9
        
        # 根据过程类型计算合频/倍频波长
        if process_type == 'SHG':
            # 二次谐波：λ₃ = λ₁ / 2
            self.lam3 = self.lam1 / 2
        else:
            # 和频：1/λ₃ = 1/λ₁ + 1/λ₂
            self.lam3 = 1 / (1/self.lam1 + 1/self.lam2)
        
        self.omega1 = 2 * np.pi * self.c / self.lam1
        self.omega2 = 2 * np.pi * self.c / self.lam2
        self.omega3 = 2 * np.pi * self.c / self.lam3
        
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3

    def _power_to_amplitude(self, P, n, w0_mm):
        """
        将功率 P 转换为电场振幅 A (Boyd定义)
        I = 2 * n * epsilon0 * c * |A|^2
        P = I * Area
        """
        if P == 0: return 0j
        w0 = w0_mm * 1e-3
        area = np.pi * w0**2
        intensity = P / area
        # |A| = sqrt( I / (2 * n * eps0 * c) )
        A_abs = np.sqrt(intensity / (2 * n * self.epsilon0 * self.c))
        return A_abs + 0j # 返回复数形式

    def _intensity_to_amplitude(self, I, n):
        """
        将强度 I 转换为电场振幅 A (Boyd定义)
        I = 2 * n * epsilon0 * c * |A|^2
        """
        if I <= 0: return 0j
        A_abs = np.sqrt(I / (2 * n * self.epsilon0 * self.c))
        return A_abs + 0j

    def _amplitude_to_intensity(self, A, n):
        """将电场振幅 A 转换回强度 I"""
        return 2 * n * self.epsilon0 * self.c * np.abs(A)**2

    def _amplitude_to_power(self, A, n, w0_mm):
        """将电场振幅 A 转换回 功率 P"""
        w0 = w0_mm * 1e-3
        area = np.pi * w0**2
        I = 2 * n * self.epsilon0 * self.c * np.abs(A)**2
        return I * area

    def solve_at_intensity(self, I1, I2, I3, delta_k, alphas):
        """
        在给定强度下求解耦合波方程
        :param I1, I2, I3: 输入强度 (W/m^2)
        :param delta_k: 相位失配 (m^-1)
        :param alphas: 吸收系数列表 (m^-1)
        :return: 出口处的强度 I1_out, I2_out, I3_out
        """
        # 转换为振幅
        A1_0 = self._intensity_to_amplitude(I1, self.n1)
        A2_0 = self._intensity_to_amplitude(I2, self.n2)
        A3_0 = self._intensity_to_amplitude(I3, self.n3)
        y0 = [A1_0, A2_0, A3_0]
        
        # 求解 ODE (简化版，只需要出口结果)
        z_span = (0, self.L)
        
        sol = solve_ivp(
            fun=lambda z, y: self.coupled_odes(z, y, delta_k, alphas),
            t_span=z_span,
            y0=y0,
            method='RK45'
        )
        
        # 取晶体出口处的振幅并转换为强度
        A1_out = sol.y[0, -1]
        A2_out = sol.y[1, -1]
        A3_out = sol.y[2, -1]
        
        I1_out = self._amplitude_to_intensity(A1_out, self.n1)
        I2_out = self._amplitude_to_intensity(A2_out, self.n2)
        I3_out = self._amplitude_to_intensity(A3_out, self.n3)
        
        return I1_out, I2_out, I3_out

    def coupled_odes(self, z, A, delta_k, alphas):
        """
        耦合波方程组 (The System of ODEs)
        输入向量 A = [A1, A2, A3] (复数)
        """
        A1, A2, A3 = A
        alpha1, alpha2, alpha3 = alphas
        
        # 根据过程类型确定耦合系数
        if self.process_type == 'SHG':
            # SHG: 没有系数2，K = i*omega*deff/(n*c)
            kappa = self.deff / self.c
        else:
            # SFG: 有系数2，K = i*2*omega*deff/(n*c)
            kappa = 2 * self.deff / self.c
        
        # 定义耦合系数 K1, K2, K3
        K1 = 1j * (self.omega1 * kappa) / self.n1
        K2 = 1j * (self.omega2 * kappa) / self.n2
        K3 = 1j * (self.omega3 * kappa) / self.n3
        
        # 相位因子
        exp_minus = np.exp(-1j * delta_k * z)
        exp_plus  = np.exp(+1j * delta_k * z)
        
        # --- 核心方程 ---
        # 注意吸收项是 alpha/2
        dA1dz = K1 * A3 * np.conj(A2) * exp_minus - (alpha1 / 2) * A1
        dA2dz = K2 * A3 * np.conj(A1) * exp_minus - (alpha2 / 2) * A2
        dA3dz = K3 * A1 * A2          * exp_plus  - (alpha3 / 2) * A3
        
        return [dA1dz, dA2dz, dA3dz]

    def solve(self, P1_W, P2_W, P3_W, w0_mm, delta_k_inv_cm, alpha_cm_inv=[0,0,0]):
        """
        执行仿真
        :param P1_W, P2_W, P3_W: 输入功率 (Watts)
        :param w0_mm: 光斑半径 (mm)
        :param delta_k_inv_cm: 相位失配 (cm^-1)
        :param alpha_cm_inv: [alpha1, alpha2, alpha3] 吸收系数列表 (cm^-1)
        """
        # 1. 单位转换
        dk = delta_k_inv_cm * 100  # cm^-1 -> m^-1
        # 吸收系数: cm^-1 -> m^-1
        alphas = [a * 100 for a in alpha_cm_inv] 
        
        # 2. 初始条件 (Power -> Amplitude)
        A1_0 = self._power_to_amplitude(P1_W, self.n1, w0_mm)
        A2_0 = self._power_to_amplitude(P2_W, self.n2, w0_mm)
        A3_0 = self._power_to_amplitude(P3_W, self.n3, w0_mm)
        y0 = [A1_0, A2_0, A3_0]
        
        # 3. 求解 ODE
        z_span = (0, self.L)
        z_eval = np.linspace(0, self.L, 500) # 仿真步数 500
        
        sol = solve_ivp(
            fun=lambda z, y: self.coupled_odes(z, y, dk, alphas),
            t_span=z_span,
            y0=y0,
            t_eval=z_eval,
            method='RK45' # Runge-Kutta 4(5)
        )
        
        # 4. 结果处理 (Amplitude -> Power)
        # sol.y 的形状是 (3, N)
        P1_evol = self._amplitude_to_power(sol.y[0], self.n1, w0_mm)
        P2_evol = self._amplitude_to_power(sol.y[1], self.n2, w0_mm)
        P3_evol = self._amplitude_to_power(sol.y[2], self.n3, w0_mm)
        
        return z_eval, P1_evol, P2_evol, P3_evol

    def solve_spatial(self, P1_W, P2_W, P3_W, w0_mm, delta_k_inv_cm, alpha_cm_inv=[0,0,0],
                      beam_profile='flat', n_radial_slices=30, n_z_points=100):
        """
        考虑空间光斑结构的求解 (空间积分)
        
        :param P1_W, P2_W, P3_W: 输入总功率 (Watts)
        :param w0_mm: 光斑半径 (mm) - 对于平顶是半径，对于高斯是1/e^2半径
        :param delta_k_inv_cm: 相位失配 (cm^-1)
        :param alpha_cm_inv: 吸收系数列表 (cm^-1)
        :param beam_profile: 光斑类型 ('flat' 平顶 或 'gaussian' 高斯)
        :param n_radial_slices: 径向切片数量 (用于空间积分)
        :param n_z_points: z方向点数 (用于热力图)
        
        :return: 包含空间分辨结果的字典
        """
        # 单位转换
        dk = delta_k_inv_cm * 100  # cm^-1 -> m^-1
        alphas = [a * 100 for a in alpha_cm_inv]  # cm^-1 -> m^-1
        w0 = w0_mm * 1e-3  # mm -> m
        
        # z轴网格
        z_eval = np.linspace(0, self.L, n_z_points)
        
        if beam_profile == 'flat':
            # 平顶光斑：整个光斑功率密度均匀
            # 计算均匀强度
            I1_uniform = P1_W / (np.pi * w0**2) if P1_W > 0 else 0
            I2_uniform = P2_W / (np.pi * w0**2) if P2_W > 0 else 0
            I3_uniform = P3_W / (np.pi * w0**2) if P3_W > 0 else 0
            
            # 求解中心点的完整演化
            A1_0 = self._intensity_to_amplitude(I1_uniform, self.n1)
            A2_0 = self._intensity_to_amplitude(I2_uniform, self.n2)
            A3_0 = self._intensity_to_amplitude(I3_uniform, self.n3)
            y0 = [A1_0, A2_0, A3_0]
            
            sol = solve_ivp(
                fun=lambda z, y: self.coupled_odes(z, y, dk, alphas),
                t_span=(0, self.L),
                y0=y0,
                t_eval=z_eval,
                method='RK45'
            )
            
            I1_z = self._amplitude_to_intensity(sol.y[0], self.n1)
            I2_z = self._amplitude_to_intensity(sol.y[1], self.n2)
            I3_z = self._amplitude_to_intensity(sol.y[2], self.n3)
            
            # 生成空间分布 (平顶: r < w0 内均匀)
            r_axis = np.linspace(0, w0 * 1.5, n_radial_slices)
            
            # 创建2D数组 [n_z, n_r]
            I1_2d = np.zeros((n_z_points, n_radial_slices))
            I2_2d = np.zeros((n_z_points, n_radial_slices))
            I3_2d = np.zeros((n_z_points, n_radial_slices))
            
            for iz in range(n_z_points):
                for ir in range(n_radial_slices):
                    if r_axis[ir] <= w0:
                        I1_2d[iz, ir] = I1_z[iz]
                        I2_2d[iz, ir] = I2_z[iz]
                        I3_2d[iz, ir] = I3_z[iz]
            
            P1_out = I1_z[-1] * np.pi * w0**2
            P2_out = I2_z[-1] * np.pi * w0**2
            P3_out = I3_z[-1] * np.pi * w0**2
            
            results = {
                'z_axis': z_eval,
                'r_axis': r_axis,
                'I1_2d': I1_2d,  # 2D强度分布 [z, r]
                'I2_2d': I2_2d,
                'I3_2d': I3_2d,
                'P1_evol': I1_z,  # 中心强度演化
                'P2_evol': I2_z,
                'P3_evol': I3_z,
                'P1_in': P1_W,
                'P2_in': P2_W,
                'P3_in': P3_W,
                'P1_out': P1_out,
                'P2_out': P2_out,
                'P3_out': P3_out,
                'beam_profile': beam_profile,
                'w0': w0,
            }
            return results
        
        else:  # gaussian
            # 高斯光斑：I(r) = I0 * exp(-2*r^2/w0^2)
            I1_peak = 2 * P1_W / (np.pi * w0**2) if P1_W > 0 else 0
            I2_peak = 2 * P2_W / (np.pi * w0**2) if P2_W > 0 else 0
            I3_peak = 2 * P3_W / (np.pi * w0**2) if P3_W > 0 else 0
            
            # 建立径向网格 (覆盖到2.5*w0)
            r_max = 2.5 * w0
            r_axis = np.linspace(0, r_max, n_radial_slices)
            
            # 创建2D数组存储完整的z-r演化 [n_z, n_r]
            I1_2d = np.zeros((n_z_points, n_radial_slices))
            I2_2d = np.zeros((n_z_points, n_radial_slices))
            I3_2d = np.zeros((n_z_points, n_radial_slices))
            
            # 对每个径向位置求解完整的z演化
            for ir in range(n_radial_slices):
                r = r_axis[ir]
                # 该半径处的输入强度
                I1_in = I1_peak * np.exp(-2 * r**2 / w0**2)
                I2_in = I2_peak * np.exp(-2 * r**2 / w0**2)
                I3_in = I3_peak * np.exp(-2 * r**2 / w0**2)
                
                if I1_in > 1e-10 or I2_in > 1e-10:
                    # 求解该径向位置的完整z演化
                    A1_0 = self._intensity_to_amplitude(I1_in, self.n1)
                    A2_0 = self._intensity_to_amplitude(I2_in, self.n2)
                    A3_0 = self._intensity_to_amplitude(I3_in, self.n3)
                    y0 = [A1_0, A2_0, A3_0]
                    
                    sol = solve_ivp(
                        fun=lambda z, y: self.coupled_odes(z, y, dk, alphas),
                        t_span=(0, self.L),
                        y0=y0,
                        t_eval=z_eval,
                        method='RK45'
                    )
                    
                    I1_2d[:, ir] = self._amplitude_to_intensity(sol.y[0], self.n1)
                    I2_2d[:, ir] = self._amplitude_to_intensity(sol.y[1], self.n2)
                    I3_2d[:, ir] = self._amplitude_to_intensity(sol.y[2], self.n3)
            
            # 计算各z位置的总功率 (空间积分)
            P1_z = np.zeros(n_z_points)
            P2_z = np.zeros(n_z_points)
            P3_z = np.zeros(n_z_points)
            
            for iz in range(n_z_points):
                P1_z[iz] = 2 * np.pi * np.trapz(I1_2d[iz, :] * r_axis, r_axis)
                P2_z[iz] = 2 * np.pi * np.trapz(I2_2d[iz, :] * r_axis, r_axis)
                P3_z[iz] = 2 * np.pi * np.trapz(I3_2d[iz, :] * r_axis, r_axis)
            
            results = {
                'z_axis': z_eval,
                'r_axis': r_axis,
                'I1_2d': I1_2d,  # 2D强度分布 [z, r]
                'I2_2d': I2_2d,
                'I3_2d': I3_2d,
                'P1_evol': I1_2d[:, 0],  # 中心强度演化 (r=0)
                'P2_evol': I2_2d[:, 0],
                'P3_evol': I3_2d[:, 0],
                'P1_z': P1_z,  # 各z位置的总功率
                'P2_z': P2_z,
                'P3_z': P3_z,
                'P1_in': P1_W,
                'P2_in': P2_W,
                'P3_in': P3_W,
                'P1_out': P1_z[-1],
                'P2_out': P2_z[-1],
                'P3_out': P3_z[-1],
                'beam_profile': beam_profile,
                'w0': w0,
                'I1_peak_in': I1_peak,
                'I2_peak_in': I2_peak,
                'I3_peak_in': I3_peak,
            }
            return results

    def plot_beam_heatmap(self, results, z_position_mm, wave='P3'):
        """
        绘制指定z位置的光斑热力图
        
        :param results: solve_spatial返回的结果字典
        :param z_position_mm: 观察位置 (mm)
        :param wave: 要显示的波 ('P1', 'P2', 'P3')
        :return: matplotlib figure
        """
        z_axis = results['z_axis']
        r_axis = results['r_axis']
        w0 = results['w0']
        
        # 找到最接近的z索引
        z_position = z_position_mm * 1e-3  # mm -> m
        iz = np.argmin(np.abs(z_axis - z_position))
        actual_z_mm = z_axis[iz] * 1e3
        
        # 获取该z位置的径向强度分布
        if wave == 'P1':
            I_r = results['I1_2d'][iz, :]
            wavelength = self.lam1 * 1e9
            label = 'Fundamental' if self.process_type == 'SHG' else 'Signal (P1)'
        elif wave == 'P2':
            I_r = results['I2_2d'][iz, :]
            wavelength = self.lam2 * 1e9
            label = 'Idler (P2)'
        else:  # P3
            I_r = results['I3_2d'][iz, :]
            wavelength = self.lam3 * 1e9
            label = 'SH' if self.process_type == 'SHG' else 'Sum (P3)'
        
        # 使用jet颜色映射
        cmap = 'jet'
        
        # 创建2D圆形热力图
        # 生成x-y网格
        r_mm = r_axis * 1e3
        n_points = 200
        x = np.linspace(-r_mm[-1], r_mm[-1], n_points)
        y = np.linspace(-r_mm[-1], r_mm[-1], n_points)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        # 插值得到2D强度分布
        from scipy.interpolate import interp1d
        interp_func = interp1d(r_mm, I_r, kind='linear', bounds_error=False, fill_value=0)
        I_2d = interp_func(R)
        
        # 转换为 W/cm²
        I_2d_cm2 = I_2d / 1e4
        
        fig, ax = plt.subplots(figsize=(8, 7))
        
        # 绘制热力图
        im = ax.pcolormesh(X, Y, I_2d_cm2, cmap=cmap, shading='auto')
        cbar = plt.colorbar(im, ax=ax, label='Intensity (W/cm²)')
        
        # 添加1/e²半径圆圈参考
        w0_mm = w0 * 1e3
        circle = plt.Circle((0, 0), w0_mm, fill=False, color='white', 
                            linestyle='--', linewidth=1.5, label=f'w₀={w0_mm:.3f}mm')
        ax.add_patch(circle)
        
        ax.set_xlabel('x (mm)', fontsize=12)
        ax.set_ylabel('y (mm)', fontsize=12)
        ax.set_title(f'{label} ({wavelength:.1f}nm) @ z={actual_z_mm:.2f}mm\n{results["beam_profile"]} beam', fontsize=14)
        ax.set_aspect('equal')
        ax.legend(loc='upper right')
        
        # 添加峰值强度标注
        peak_I = np.max(I_2d_cm2)
        ax.text(0.02, 0.98, f'Peak: {peak_I:.2e} W/cm²', transform=ax.transAxes, 
                fontsize=10, verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        fig.tight_layout()
        return fig

    def plot_results(self, z, P1, P2, P3):
        fig, ax = plt.subplots(figsize=(10, 6))
        # 转换 z 轴为 mm 方便查看
        z_mm = z * 1e3
        
        if self.process_type == 'SHG':
            # SHG只显示基频和倍频
            ax.plot(z_mm, P1, 'b-', label=f'Fundamental (P1): {self.lam1*1e9:.1f}nm', linewidth=2)
            ax.plot(z_mm, P3, 'g-', label=f'SH (P3): {self.lam3*1e9:.1f}nm', linewidth=2)
        else:
            # SFG显示三束光
            ax.plot(z_mm, P1, label=f'Signal (P1): {self.lam1*1e9:.1f}nm', linewidth=2)
            ax.plot(z_mm, P2, label=f'Idler (P2): {self.lam2*1e9:.1f}nm', linewidth=2)
            ax.plot(z_mm, P3, label=f'Sum (P3): {self.lam3*1e9:.1f}nm', linewidth=2)
        
        ax.set_xlabel('Crystal Length (mm)', fontsize=12)
        ax.set_ylabel('Power (W)', fontsize=12)
        title = f'{self.process_type} Power Evolution\n(deff={self.deff*1e12:.1f} pm/V)'
        ax.set_title(title, fontsize=14)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()
        return fig

    def solve_pulse(self, P1_avg, P2_avg, P3_avg, w0_mm, delta_k_inv_cm, 
                    alpha_cm_inv=[0,0,0], rep_rate_Hz=1e4, pulse_width_s=10e-9, 
                    n_time_slices=100, pulse_shape='gaussian'):
        """
        脉冲光的切片法求解 (Quasi-Static Approximation)
        
        :param P1_avg, P2_avg, P3_avg: 平均功率 (W)
        :param w0_mm: 光斑半径 (mm)
        :param delta_k_inv_cm: 相位失配 (cm^-1)
        :param alpha_cm_inv: 吸收系数列表 (cm^-1)
        :param rep_rate_Hz: 重复频率 (Hz)
        :param pulse_width_s: 脉冲宽度 FWHM (s)
        :param n_time_slices: 时间切片数量
        :param pulse_shape: 脉冲形状 ('gaussian' 或 'square')
        
        :return: t_axis, P1_in, P2_in, P3_in, P1_out, P2_out, P3_out, pulse_energy_in, pulse_energy_out
        """
        # 1. 计算脉冲能量和峰值功率
        duty_cycle = rep_rate_Hz * pulse_width_s
        
        # 单脉冲能量 = 平均功率 / 重频
        E1_pulse = P1_avg / rep_rate_Hz if rep_rate_Hz > 0 else 0  # J
        E2_pulse = P2_avg / rep_rate_Hz if rep_rate_Hz > 0 else 0
        E3_pulse = P3_avg / rep_rate_Hz if rep_rate_Hz > 0 else 0
        
        # 2. 建立时间轴 (以脉冲中心为原点，覆盖 ±3σ)
        if pulse_shape == 'gaussian':
            # 高斯脉冲: FWHM = 2.355 * sigma
            sigma = pulse_width_s / 2.355
            t_range = 6 * sigma  # 覆盖 ±3σ
        else:
            # 方波脉冲
            t_range = pulse_width_s * 1.5
        
        t_axis = np.linspace(-t_range/2, t_range/2, n_time_slices)
        dt = t_axis[1] - t_axis[0]
        
        # 3. 生成脉冲形状 (归一化使得积分=1)
        if pulse_shape == 'gaussian':
            # 高斯脉冲形状
            pulse_envelope = np.exp(-t_axis**2 / (2 * sigma**2))
            # 归一化: ∫ envelope * dt = 1
            pulse_envelope = pulse_envelope / (np.sum(pulse_envelope) * dt)
        else:
            # 方波脉冲
            pulse_envelope = np.where(np.abs(t_axis) <= pulse_width_s/2, 1.0, 0.0)
            pulse_envelope = pulse_envelope / (np.sum(pulse_envelope) * dt)
        
        # 4. 计算每个时间片的瞬时功率
        # P(t) = E_pulse * envelope(t), 其中 ∫envelope*dt = 1
        P1_in = E1_pulse * pulse_envelope
        P2_in = E2_pulse * pulse_envelope
        P3_in = E3_pulse * pulse_envelope
        
        # 5. 对每个时间片求解耦合波方程
        P1_out = np.zeros_like(t_axis)
        P2_out = np.zeros_like(t_axis)
        P3_out = np.zeros_like(t_axis)
        
        for i in range(n_time_slices):
            # 只对有功率的时间片求解
            if P1_in[i] > 1e-10 or P2_in[i] > 1e-10:
                z, p1, p2, p3 = self.solve(
                    P1_W=P1_in[i],
                    P2_W=P2_in[i],
                    P3_W=P3_in[i],
                    w0_mm=w0_mm,
                    delta_k_inv_cm=delta_k_inv_cm,
                    alpha_cm_inv=alpha_cm_inv
                )
                # 取晶体出口处的功率
                P1_out[i] = p1[-1]
                P2_out[i] = p2[-1]
                P3_out[i] = p3[-1]
            else:
                P1_out[i] = P1_in[i]
                P2_out[i] = P2_in[i]
                P3_out[i] = P3_in[i]
        
        # 6. 计算输出脉冲能量
        E1_out = np.trapz(P1_out, t_axis)
        E2_out = np.trapz(P2_out, t_axis)
        E3_out = np.trapz(P3_out, t_axis)
        
        # 7. 计算输出平均功率
        P1_avg_out = E1_out * rep_rate_Hz
        P2_avg_out = E2_out * rep_rate_Hz
        P3_avg_out = E3_out * rep_rate_Hz
        
        results = {
            't_axis': t_axis,           # 时间轴 (s)
            'P1_in': P1_in,             # 输入脉冲形状 P1(t)
            'P2_in': P2_in,
            'P3_in': P3_in,
            'P1_out': P1_out,           # 输出脉冲形状 P1(t)
            'P2_out': P2_out,
            'P3_out': P3_out,
            'E1_in': E1_pulse,          # 输入单脉冲能量 (J)
            'E2_in': E2_pulse,
            'E3_in': E3_pulse,
            'E1_out': E1_out,           # 输出单脉冲能量 (J)
            'E2_out': E2_out,
            'E3_out': E3_out,
            'P1_avg_out': P1_avg_out,   # 输出平均功率 (W)
            'P2_avg_out': P2_avg_out,
            'P3_avg_out': P3_avg_out,
            'P1_peak_in': np.max(P1_in),    # 输入峰值功率 (W)
            'P2_peak_in': np.max(P2_in),
            'P3_peak_in': np.max(P3_in),
            'P1_peak_out': np.max(P1_out),  # 输出峰值功率 (W)
            'P2_peak_out': np.max(P2_out),
            'P3_peak_out': np.max(P3_out),
        }
        
        return results

    def plot_pulse_results(self, results):
        """绘制脉冲仿真结果"""
        t_ns = results['t_axis'] * 1e9  # 转换为 ns
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # 左图：输入脉冲形状
        ax1 = axes[0]
        if self.process_type == 'SHG':
            ax1.plot(t_ns, results['P1_in']/1e3, 'b-', label=f'Fundamental: {self.lam1*1e9:.1f}nm', linewidth=2)
            ax1.plot(t_ns, results['P3_in']/1e3, 'g-', label=f'SH: {self.lam3*1e9:.1f}nm', linewidth=2)
        else:
            ax1.plot(t_ns, results['P1_in']/1e3, 'b-', label=f'P1: {self.lam1*1e9:.1f}nm', linewidth=2)
            ax1.plot(t_ns, results['P2_in']/1e3, 'orange', label=f'P2: {self.lam2*1e9:.1f}nm', linewidth=2)
            ax1.plot(t_ns, results['P3_in']/1e3, 'g-', label=f'P3: {self.lam3*1e9:.1f}nm', linewidth=2)
        ax1.set_xlabel('Time (ns)', fontsize=12)
        ax1.set_ylabel('Instantaneous Power (kW)', fontsize=12)
        ax1.set_title('Input Pulse Shape', fontsize=14)
        ax1.legend()
        ax1.grid(True, linestyle='--', alpha=0.6)
        
        # 右图：输出脉冲形状
        ax2 = axes[1]
        if self.process_type == 'SHG':
            ax2.plot(t_ns, results['P1_out']/1e3, 'b-', label=f'Fundamental: {self.lam1*1e9:.1f}nm', linewidth=2)
            ax2.plot(t_ns, results['P3_out']/1e3, 'g-', label=f'SH: {self.lam3*1e9:.1f}nm', linewidth=2)
        else:
            ax2.plot(t_ns, results['P1_out']/1e3, 'b-', label=f'P1: {self.lam1*1e9:.1f}nm', linewidth=2)
            ax2.plot(t_ns, results['P2_out']/1e3, 'orange', label=f'P2: {self.lam2*1e9:.1f}nm', linewidth=2)
            ax2.plot(t_ns, results['P3_out']/1e3, 'g-', label=f'P3: {self.lam3*1e9:.1f}nm', linewidth=2)
        ax2.set_xlabel('Time (ns)', fontsize=12)
        ax2.set_ylabel('Instantaneous Power (kW)', fontsize=12)
        ax2.set_title(f'Output Pulse Shape ({self.process_type})', fontsize=14)
        ax2.legend()
        ax2.grid(True, linestyle='--', alpha=0.6)
        
        fig.tight_layout()
        return fig

    def solve_pulse_spatial(self, P1_avg, P2_avg, P3_avg, w0_mm, delta_k_inv_cm, 
                            alpha_cm_inv=[0,0,0], rep_rate_Hz=1e4, pulse_width_s=10e-9, 
                            n_time_slices=50, n_radial_slices=30, pulse_shape='gaussian',
                            beam_profile='gaussian'):
        """
        脉冲光的时间-空间联合求解
        
        对于每个时间切片，计算空间（径向）分布，然后时间积分得到能量密度分布
        
        :param P1_avg, P2_avg, P3_avg: 平均功率 (W)
        :param w0_mm: 光斑半径 (mm)
        :param delta_k_inv_cm: 相位失配 (cm^-1)
        :param alpha_cm_inv: 吸收系数列表 (cm^-1)
        :param rep_rate_Hz: 重复频率 (Hz)
        :param pulse_width_s: 脉冲宽度 FWHM (s)
        :param n_time_slices: 时间切片数量
        :param n_radial_slices: 径向切片数量
        :param pulse_shape: 脉冲时间形状 ('gaussian' 或 'square')
        :param beam_profile: 空间光斑形状 ('gaussian' 或 'flat-top')
        
        :return: 包含时间-空间演化结果的字典
        """
        w0 = w0_mm * 1e-3  # mm -> m
        dk = delta_k_inv_cm * 1e2  # cm^-1 -> m^-1
        alphas = [a * 1e2 for a in alpha_cm_inv]  # cm^-1 -> m^-1
        
        # 1. 计算单脉冲能量
        E1_pulse = P1_avg / rep_rate_Hz if rep_rate_Hz > 0 else 0  # J
        E2_pulse = P2_avg / rep_rate_Hz if rep_rate_Hz > 0 else 0
        E3_pulse = P3_avg / rep_rate_Hz if rep_rate_Hz > 0 else 0
        
        # 2. 建立时间轴
        if pulse_shape == 'gaussian':
            sigma_t = pulse_width_s / 2.355
            t_range = 6 * sigma_t
        else:
            t_range = pulse_width_s * 1.5
        
        t_axis = np.linspace(-t_range/2, t_range/2, n_time_slices)
        dt = t_axis[1] - t_axis[0]
        
        # 3. 生成脉冲时间包络 (归一化使得 ∫envelope*dt = 1)
        if pulse_shape == 'gaussian':
            pulse_envelope = np.exp(-t_axis**2 / (2 * sigma_t**2))
            pulse_envelope = pulse_envelope / (np.sum(pulse_envelope) * dt)
        else:
            pulse_envelope = np.where(np.abs(t_axis) <= pulse_width_s/2, 1.0, 0.0)
            pulse_envelope = pulse_envelope / (np.sum(pulse_envelope) * dt)
        
        # 4. 建立径向网格
        if beam_profile == 'flat-top':
            r_max = w0 * 1.5
        else:
            r_max = w0 * 2.5
        r_axis = np.linspace(0, r_max, n_radial_slices)
        
        # 5. 建立z轴
        n_z_points = 100
        z_axis = np.linspace(0, self.L, n_z_points)
        
        # 6. 创建3D数组 [t, z, r] 存储瞬时强度
        I1_3d = np.zeros((n_time_slices, n_z_points, n_radial_slices))
        I2_3d = np.zeros((n_time_slices, n_z_points, n_radial_slices))
        I3_3d = np.zeros((n_time_slices, n_z_points, n_radial_slices))
        
        # 7. 计算每个时间切片的瞬时功率
        P1_t = E1_pulse * pulse_envelope  # W, 瞬时功率随时间分布
        P2_t = E2_pulse * pulse_envelope
        P3_t = E3_pulse * pulse_envelope
        
        # 8. 对每个时间切片求解空间演化
        for it in range(n_time_slices):
            P1_inst = P1_t[it]
            P2_inst = P2_t[it]
            P3_inst = P3_t[it]
            
            if P1_inst < 1e-10 and P2_inst < 1e-10:
                continue
            
            if beam_profile == 'flat-top':
                # 平顶光斑：均匀强度
                I1_in = P1_inst / (np.pi * w0**2)
                I2_in = P2_inst / (np.pi * w0**2)
                I3_in = P3_inst / (np.pi * w0**2)
                
                # 求解z演化
                A1_0 = self._intensity_to_amplitude(I1_in, self.n1)
                A2_0 = self._intensity_to_amplitude(I2_in, self.n2)
                A3_0 = self._intensity_to_amplitude(I3_in, self.n3)
                y0 = [A1_0, A2_0, A3_0]
                
                sol = solve_ivp(
                    fun=lambda z, y: self.coupled_odes(z, y, dk, alphas),
                    t_span=(0, self.L),
                    y0=y0,
                    t_eval=z_axis,
                    method='RK45'
                )
                
                I1_z = self._amplitude_to_intensity(sol.y[0], self.n1)
                I2_z = self._amplitude_to_intensity(sol.y[1], self.n2)
                I3_z = self._amplitude_to_intensity(sol.y[2], self.n3)
                
                # 填充该时间切片的空间分布
                for ir in range(n_radial_slices):
                    if r_axis[ir] <= w0:
                        I1_3d[it, :, ir] = I1_z
                        I2_3d[it, :, ir] = I2_z
                        I3_3d[it, :, ir] = I3_z
            
            else:  # gaussian beam
                # 高斯光斑：I(r) = I0 * exp(-2*r^2/w0^2)
                I1_peak = 2 * P1_inst / (np.pi * w0**2)
                I2_peak = 2 * P2_inst / (np.pi * w0**2)
                I3_peak = 2 * P3_inst / (np.pi * w0**2)
                
                # 对每个径向位置求解
                for ir in range(n_radial_slices):
                    r = r_axis[ir]
                    I1_in = I1_peak * np.exp(-2 * r**2 / w0**2)
                    I2_in = I2_peak * np.exp(-2 * r**2 / w0**2)
                    I3_in = I3_peak * np.exp(-2 * r**2 / w0**2)
                    
                    if I1_in > 1e-10 or I2_in > 1e-10:
                        A1_0 = self._intensity_to_amplitude(I1_in, self.n1)
                        A2_0 = self._intensity_to_amplitude(I2_in, self.n2)
                        A3_0 = self._intensity_to_amplitude(I3_in, self.n3)
                        y0 = [A1_0, A2_0, A3_0]
                        
                        sol = solve_ivp(
                            fun=lambda z, y: self.coupled_odes(z, y, dk, alphas),
                            t_span=(0, self.L),
                            y0=y0,
                            t_eval=z_axis,
                            method='RK45'
                        )
                        
                        I1_3d[it, :, ir] = self._amplitude_to_intensity(sol.y[0], self.n1)
                        I2_3d[it, :, ir] = self._amplitude_to_intensity(sol.y[1], self.n2)
                        I3_3d[it, :, ir] = self._amplitude_to_intensity(sol.y[2], self.n3)
        
        # 9. 时间积分得到能量密度分布 F(z,r) = ∫I(t,z,r)*dt [J/m²]
        F1_2d = np.trapz(I1_3d, t_axis, axis=0)  # [z, r]
        F2_2d = np.trapz(I2_3d, t_axis, axis=0)
        F3_2d = np.trapz(I3_3d, t_axis, axis=0)
        
        # 10. 空间积分得到各z位置的脉冲能量 E(z) = 2π∫F(z,r)*r*dr [J]
        E1_z = np.zeros(n_z_points)
        E2_z = np.zeros(n_z_points)
        E3_z = np.zeros(n_z_points)
        
        for iz in range(n_z_points):
            E1_z[iz] = 2 * np.pi * np.trapz(F1_2d[iz, :] * r_axis, r_axis)
            E2_z[iz] = 2 * np.pi * np.trapz(F2_2d[iz, :] * r_axis, r_axis)
            E3_z[iz] = 2 * np.pi * np.trapz(F3_2d[iz, :] * r_axis, r_axis)
        
        # 11. 输出结果
        results = {
            't_axis': t_axis,           # 时间轴 (s)
            'z_axis': z_axis,           # 位置轴 (m)
            'r_axis': r_axis,           # 径向轴 (m)
            'w0': w0,                   # 光斑半径 (m)
            'beam_profile': beam_profile,
            'pulse_shape': pulse_shape,
            
            # 3D瞬时强度分布 [t, z, r]
            'I1_3d': I1_3d,
            'I2_3d': I2_3d,
            'I3_3d': I3_3d,
            
            # 输入/输出能量密度分布 [r] (J/m²)
            'F1_in': F1_2d[0, :],       # 输入能量密度径向分布
            'F2_in': F2_2d[0, :],
            'F3_in': F3_2d[0, :],
            'F1_out': F1_2d[-1, :],     # 输出能量密度径向分布
            'F2_out': F2_2d[-1, :],
            'F3_out': F3_2d[-1, :],
            
            # 2D能量密度分布 [z, r] (J/m²)
            'F1_2d': F1_2d,
            'F2_2d': F2_2d,
            'F3_2d': F3_2d,
            
            # 各z位置的脉冲能量 (J)
            'E1_z': E1_z,
            'E2_z': E2_z,
            'E3_z': E3_z,
            
            # 输入/输出脉冲能量 (J)
            'E1_in': E1_pulse,
            'E2_in': E2_pulse,
            'E3_in': E3_pulse,
            'E1_out': E1_z[-1],
            'E2_out': E2_z[-1],
            'E3_out': E3_z[-1],
            
            # 时间切片信息
            'P1_t': P1_t,               # 各时间点的瞬时功率 (W)
            'P2_t': P2_t,
            'P3_t': P3_t,
            'P1_peak_in': np.max(P1_t),
            'P2_peak_in': np.max(P2_t),
            'P3_peak_in': np.max(P3_t),
        }
        
        return results

    def plot_pulse_beam_heatmap(self, results, wave='P3', position='output'):
        """
        绘制脉冲输出光斑的能量密度热力图
        
        :param results: solve_pulse_spatial返回的结果字典
        :param wave: 要显示的波 ('P1', 'P2', 'P3')
        :param position: 'input' 或 'output'
        :return: matplotlib figure
        """
        r_axis = results['r_axis']
        w0 = results['w0']
        beam_profile = results['beam_profile']
        
        # 获取能量密度分布
        if wave == 'P1':
            F_r = results['F1_out'] if position == 'output' else results['F1_in']
            wavelength = self.lam1 * 1e9
            label = 'Fundamental' if self.process_type == 'SHG' else 'Signal (P1)'
        elif wave == 'P2':
            F_r = results['F2_out'] if position == 'output' else results['F2_in']
            wavelength = self.lam2 * 1e9
            label = 'Idler (P2)'
        else:  # P3
            F_r = results['F3_out'] if position == 'output' else results['F3_in']
            wavelength = self.lam3 * 1e9
            label = 'SH' if self.process_type == 'SHG' else 'Sum (P3)'
        
        # 创建2D圆形热力图
        r_mm = r_axis * 1e3
        n_points = 200
        x = np.linspace(-r_mm[-1], r_mm[-1], n_points)
        y = np.linspace(-r_mm[-1], r_mm[-1], n_points)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        # 插值
        from scipy.interpolate import interp1d
        interp_func = interp1d(r_mm, F_r, kind='linear', bounds_error=False, fill_value=0)
        F_2d = interp_func(R)
        
        # 转换为 J/cm² (常用单位)
        F_2d_cm2 = F_2d / 1e4
        
        fig, ax = plt.subplots(figsize=(8, 7))
        
        im = ax.pcolormesh(X, Y, F_2d_cm2, cmap='jet', shading='auto')
        cbar = plt.colorbar(im, ax=ax, label='Fluence (J/cm²)')
        
        # 添加1/e²半径圆圈
        w0_mm = w0 * 1e3
        circle = plt.Circle((0, 0), w0_mm, fill=False, color='white', 
                            linestyle='--', linewidth=1.5, label=f'w₀={w0_mm:.3f}mm')
        ax.add_patch(circle)
        
        pos_label = 'Output' if position == 'output' else 'Input'
        ax.set_xlabel('x (mm)', fontsize=12)
        ax.set_ylabel('y (mm)', fontsize=12)
        ax.set_title(f'{pos_label} {label} ({wavelength:.1f}nm)\nPulse Fluence ({beam_profile} beam)', fontsize=14)
        ax.set_aspect('equal')
        ax.legend(loc='upper right')
        
        # 添加峰值能量密度标注
        peak_F = np.max(F_2d_cm2)
        ax.text(0.02, 0.98, f'Peak: {peak_F:.2e} J/cm²', transform=ax.transAxes, 
                fontsize=10, verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        fig.tight_layout()
        return fig

    def plot_pulse_beam_comparison(self, results, wave='P3'):
        """
        并排对比输入/输出光斑的能量密度分布
        
        :param results: solve_pulse_spatial返回的结果字典
        :param wave: 要显示的波 ('P1', 'P2', 'P3')
        :return: matplotlib figure
        """
        r_axis = results['r_axis']
        w0 = results['w0']
        
        if wave == 'P1':
            F_in = results['F1_in']
            F_out = results['F1_out']
            wavelength = self.lam1 * 1e9
            label = 'Fundamental' if self.process_type == 'SHG' else 'Signal (P1)'
        elif wave == 'P2':
            F_in = results['F2_in']
            F_out = results['F2_out']
            wavelength = self.lam2 * 1e9
            label = 'Idler (P2)'
        else:
            F_in = results['F3_in']
            F_out = results['F3_out']
            wavelength = self.lam3 * 1e9
            label = 'SH' if self.process_type == 'SHG' else 'Sum (P3)'
        
        # 创建2D网格
        r_mm = r_axis * 1e3
        n_points = 200
        x = np.linspace(-r_mm[-1], r_mm[-1], n_points)
        y = np.linspace(-r_mm[-1], r_mm[-1], n_points)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        from scipy.interpolate import interp1d
        interp_in = interp1d(r_mm, F_in, kind='linear', bounds_error=False, fill_value=0)
        interp_out = interp1d(r_mm, F_out, kind='linear', bounds_error=False, fill_value=0)
        
        F_in_2d = interp_in(R) / 1e4  # J/cm²
        F_out_2d = interp_out(R) / 1e4
        
        # 统一色标范围
        vmax = max(np.max(F_in_2d), np.max(F_out_2d))
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # 输入光斑
        ax1 = axes[0]
        im1 = ax1.pcolormesh(X, Y, F_in_2d, cmap='jet', shading='auto', vmin=0, vmax=vmax)
        plt.colorbar(im1, ax=ax1, label='Fluence (J/cm²)')
        w0_mm = w0 * 1e3
        circle1 = plt.Circle((0, 0), w0_mm, fill=False, color='white', linestyle='--', linewidth=1.5)
        ax1.add_patch(circle1)
        ax1.set_xlabel('x (mm)')
        ax1.set_ylabel('y (mm)')
        ax1.set_title(f'Input {label}\nPeak: {np.max(F_in_2d):.2e} J/cm²', fontsize=12)
        ax1.set_aspect('equal')
        
        # 输出光斑
        ax2 = axes[1]
        im2 = ax2.pcolormesh(X, Y, F_out_2d, cmap='jet', shading='auto', vmin=0, vmax=vmax)
        plt.colorbar(im2, ax=ax2, label='Fluence (J/cm²)')
        circle2 = plt.Circle((0, 0), w0_mm, fill=False, color='white', linestyle='--', linewidth=1.5)
        ax2.add_patch(circle2)
        ax2.set_xlabel('x (mm)')
        ax2.set_ylabel('y (mm)')
        ax2.set_title(f'Output {label}\nPeak: {np.max(F_out_2d):.2e} J/cm²', fontsize=12)
        ax2.set_aspect('equal')
        
        fig.suptitle(f'{label} ({wavelength:.1f}nm) Pulse Beam Profile', fontsize=14, y=1.02)
        fig.tight_layout()
        return fig

    def plot_pulse_energy_evolution(self, results):
        """
        绘制脉冲能量沿晶体的演化
        
        :param results: solve_pulse_spatial返回的结果字典
        :return: matplotlib figure
        """
        z_mm = results['z_axis'] * 1e3
        
        # 转换为 μJ
        E1_uJ = results['E1_z'] * 1e6
        E2_uJ = results['E2_z'] * 1e6
        E3_uJ = results['E3_z'] * 1e6
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        if self.process_type == 'SHG':
            ax.plot(z_mm, E1_uJ, 'b-', label=f'Fundamental (P1): {self.lam1*1e9:.1f}nm', linewidth=2)
            ax.plot(z_mm, E3_uJ, 'g-', label=f'SH (P3): {self.lam3*1e9:.1f}nm', linewidth=2)
        else:
            ax.plot(z_mm, E1_uJ, 'b-', label=f'Signal (P1): {self.lam1*1e9:.1f}nm', linewidth=2)
            ax.plot(z_mm, E2_uJ, 'orange', label=f'Idler (P2): {self.lam2*1e9:.1f}nm', linewidth=2)
            ax.plot(z_mm, E3_uJ, 'g-', label=f'Sum (P3): {self.lam3*1e9:.1f}nm', linewidth=2)
        
        ax.set_xlabel('Crystal Length (mm)', fontsize=12)
        ax.set_ylabel('Pulse Energy (μJ)', fontsize=12)
        ax.set_title(f'{self.process_type} Pulse Energy Evolution', fontsize=14)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        
        fig.tight_layout()
        return fig


