import streamlit as st
import numpy as np
from cp_wave import ThreeWaveMixing

# --- Streamlit 应用界面 ---

st.set_page_config(page_title="非线性光学耦合波方程求解器", layout="wide")
st.title(" 非线性光学耦合波方程求解器")
st.markdown("**基于Boyd体系的三波混频仿真**")

# 侧边栏参数输入
st.sidebar.header("📊 输入参数")

# 选择过程类型
st.sidebar.subheader("0️⃣ 非线性过程")
process_type = st.sidebar.radio("过程类型", ["SFG (和频)", "SHG (倍频)"], index=0)
process_type = 'SFG' if 'SFG' in process_type else 'SHG'

st.sidebar.subheader("1️⃣ 波长参数 (nm)")
if process_type == 'SHG':
    lambda1_nm = st.sidebar.number_input("λ₁ - 基频光波长", value=1064.0, min_value=100.0, max_value=10000.0, step=1.0)
    lambda2_nm = lambda1_nm  # SHG中两个输入波长相同
    st.sidebar.info(f"SHG过程：λ₂ = λ₁ = {lambda1_nm} nm")
    # 计算倍频波长
    lambda3_calc = lambda1_nm / 2
    st.sidebar.info(f"倍频波长 λ₃ = {lambda3_calc:.2f} nm")
else:
    lambda1_nm = st.sidebar.number_input("λ₁ - 信号光波长", value=1064.0, min_value=100.0, max_value=10000.0, step=1.0)
    lambda2_nm = st.sidebar.number_input("λ₂ - 闲频光波长", value=532.0, min_value=100.0, max_value=10000.0, step=1.0)
    # 计算和频波长
    lambda3_calc = 1 / (1/lambda1_nm + 1/lambda2_nm)
    st.sidebar.info(f"和频波长 λ₃ = {lambda3_calc:.2f} nm")

st.sidebar.subheader("2️⃣ 折射率")
n1 = st.sidebar.number_input("n₁", value=1.60, min_value=1.0, max_value=3.0, step=0.01, format="%.3f")
if process_type == 'SHG':
    n2 = n1  # SHG中n₂=n₁
    st.sidebar.info(f"SHG过程：n₂ = n₁ = {n1:.3f}")
else:
    n2 = st.sidebar.number_input("n₂", value=1.61, min_value=1.0, max_value=3.0, step=0.01, format="%.3f")
n3 = st.sidebar.number_input("n₃", value=1.63, min_value=1.0, max_value=3.0, step=0.01, format="%.3f")

st.sidebar.subheader("3️⃣ 晶体参数")
deff_pm_V = st.sidebar.number_input("dₑff 有效非线性系数 (pm/V)", value=0.85, min_value=0.01, max_value=100.0, step=0.01, format="%.2f")
L_mm = st.sidebar.number_input("晶体长度 (mm)", value=15.0, min_value=0.1, max_value=100.0, step=0.1, format="%.1f")

st.sidebar.subheader("4️⃣ 输入光束参数")

# 选择连续光或脉冲光
light_mode = st.sidebar.radio("光源类型", ["连续光 (CW)", "脉冲光 (Pulsed)"])

if light_mode == "连续光 (CW)":
    # 连续光：直接输入功率
    P1_W = st.sidebar.number_input("P₁ 输入功率 (W)", value=100.0, min_value=0.0, max_value=10000.0, step=1.0)
    if process_type == 'SHG':
        P2_W = P1_W  # SHG中P2=P1
    else:
        P2_W = st.sidebar.number_input("P₂ 输入功率 (W)", value=20.0, min_value=0.0, max_value=10000.0, step=1.0)
    P3_W = st.sidebar.number_input("P₃ 初始功率 (W)", value=0.0, min_value=0.0, max_value=10000.0, step=1.0)
    # 连续光的峰值功率就是输入功率
    P1_peak = P1_W
    P2_peak = P2_W
    P3_peak = P3_W
    rep_rate = None
    pulse_width = None
else:
    # 脉冲光：输入平均功率、重频、脉宽
    st.sidebar.markdown("**脉冲参数**")
    rep_rate_kHz = st.sidebar.number_input("重复频率 (kHz)", value=10.0, min_value=0.001, max_value=1000000.0, step=1.0, format="%.3f")
    pulse_width_ns = st.sidebar.number_input("脉冲宽度 FWHM (ns)", value=10.0, min_value=0.001, max_value=10000.0, step=0.1, format="%.3f")
    pulse_shape = st.sidebar.selectbox("脉冲形状", ["gaussian", "square"], index=0)
    n_time_slices = st.sidebar.slider("时间切片数", min_value=20, max_value=200, value=50, step=10)
    
    rep_rate = rep_rate_kHz * 1e3  # kHz -> Hz
    pulse_width = pulse_width_ns * 1e-9  # ns -> s
    
    st.sidebar.markdown("**平均功率**")
    P1_avg = st.sidebar.number_input("P₁ 平均功率 (W)", value=10.0, min_value=0.0, max_value=10000.0, step=0.1, format="%.3f")
    if process_type == 'SHG':
        P2_avg = P1_avg  # SHG中P2=P1
    else:
        P2_avg = st.sidebar.number_input("P₂ 平均功率 (W)", value=10.0, min_value=0.0, max_value=10000.0, step=0.1, format="%.3f")
    P3_avg = st.sidebar.number_input("P₃ 初始平均功率 (W)", value=0.0, min_value=0.0, max_value=10000.0, step=0.1, format="%.3f")
    
    # 连续光模式不使用的变量
    P1_peak = None
    P2_peak = None
    P3_peak = None

w0_mm = st.sidebar.number_input("光斑半径 w₀ (mm)", value=0.2, min_value=0.001, max_value=10.0, step=0.001, format="%.3f")

# 光斑类型选择
beam_profile = st.sidebar.selectbox("光斑类型", ["平顶光斑 (Flat-top)", "高斯光斑 (Gaussian)"], index=1)
beam_profile_code = 'flat' if '平顶' in beam_profile else 'gaussian'
if beam_profile_code == 'gaussian':
    n_radial_slices = st.sidebar.slider("径向切片数", min_value=10, max_value=100, value=30, step=5)
else:
    n_radial_slices = 30  # 平顶光斑不需要，但给个默认值

st.sidebar.subheader("5️⃣ 相位失配与吸收")
delta_k_inv_cm = st.sidebar.number_input("Δk 相位失配 (cm⁻¹)", value=0.0, min_value=-100.0, max_value=100.0, step=0.01, format="%.3f")
alpha1 = st.sidebar.number_input("α₁ 吸收系数 (cm⁻¹)", value=0.000, min_value=0.0, max_value=10.0, step=0.001, format="%.4f")
if process_type == 'SHG':
    alpha2 = 0.0  # SHG中没有P2，α₂=0
else:
    alpha2 = st.sidebar.number_input("α₂ 吸收系数 (cm⁻¹)", value=0.000, min_value=0.0, max_value=10.0, step=0.001, format="%.4f")
alpha3 = st.sidebar.number_input("α₃ 吸收系数 (cm⁻¹)", value=0.0, min_value=0.0, max_value=10.0, step=0.01, format="%.3f")

# 运行仿真按钮
if st.sidebar.button("🚀 运行仿真", type="primary"):
    try:
        with st.spinner("正在求解耦合波方程..."):
            # 初始化模拟
            sim = ThreeWaveMixing(
                lambda1_nm=lambda1_nm,
                lambda2_nm=lambda2_nm,
                n1=n1, n2=n2, n3=n3,
                deff_pm_V=deff_pm_V,
                L_mm=L_mm,
                process_type=process_type
            )
            
            # 计算输出波长
            if process_type == 'SHG':
                lam3_nm = lambda1_nm / 2
            else:
                lam3_nm = 1 / (1/lambda1_nm + 1/lambda2_nm)
            
            import matplotlib.pyplot as plt
            plt.close('all')  # 清除之前的图形缓存
            
            # 保存sim对象和参数到session_state
            st.session_state['sim'] = sim
            st.session_state['lam3_nm'] = lam3_nm
            st.session_state['L_mm'] = L_mm
            st.session_state['process_type'] = process_type
            st.session_state['beam_profile_code'] = beam_profile_code
            
            # 根据光源类型选择不同的求解方法
            if light_mode == "连续光 (CW)":
                # 连续光：使用空间积分求解
                results = sim.solve_spatial(
                    P1_W=P1_peak,
                    P2_W=P2_peak,
                    P3_W=P3_peak,
                    w0_mm=w0_mm,
                    delta_k_inv_cm=delta_k_inv_cm,
                    alpha_cm_inv=[alpha1, alpha2, alpha3],
                    beam_profile=beam_profile_code,
                    n_radial_slices=n_radial_slices
                )
                
                # 保存所有结果到session_state
                st.session_state['cw_results'] = results
                st.session_state['cw_mode'] = True
                st.session_state['P1_peak'] = P1_peak
                st.session_state['P2_peak'] = P2_peak
                st.session_state['P3_peak'] = P3_peak
                st.session_state['beam_profile'] = beam_profile
                st.session_state['run_success'] = True
            
            else:
                # 脉冲光：使用时间-空间联合求解
                results = sim.solve_pulse_spatial(
                    P1_avg=P1_avg,
                    P2_avg=P2_avg,
                    P3_avg=P3_avg,
                    w0_mm=w0_mm,
                    delta_k_inv_cm=delta_k_inv_cm,
                    alpha_cm_inv=[alpha1, alpha2, alpha3],
                    rep_rate_Hz=rep_rate,
                    pulse_width_s=pulse_width,
                    n_time_slices=n_time_slices,
                    n_radial_slices=n_radial_slices,
                    pulse_shape=pulse_shape,
                    beam_profile=beam_profile_code
                )
                
                # 保存到session_state供后续显示
                st.session_state['pulse_spatial_results'] = results
                st.session_state['pulse_mode'] = 'spatial'
                st.session_state['cw_mode'] = False
                st.session_state['pulse_beam_profile'] = beam_profile_code
                st.session_state['pulse_lam1_nm'] = sim.lam1 * 1e9
                st.session_state['pulse_lam2_nm'] = sim.lam2 * 1e9
                st.session_state['pulse_lam3_nm'] = sim.lam3 * 1e9
                st.session_state['pulse_process_type'] = process_type
                st.session_state['run_success'] = True
                
                # 显示结果
                st.success(f"✅ {process_type} 脉冲仿真完成！")
                
                process_label = "倍频" if process_type == "SHG" else "和频"
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.markdown("**波长信息**")
                    st.metric(f"{process_label}波长 λ₃", f"{lam3_nm:.2f} nm")
                    st.markdown("**脉冲参数**")
                    st.write(f"重频: {rep_rate/1e3:.3f} kHz")
                    st.write(f"脉宽 FWHM: {pulse_width*1e9:.3f} ns")
                    st.write(f"脉冲形状: {pulse_shape}")
                    st.write(f"时间切片数: {n_time_slices}")
                    st.write(f"空间光斑: {beam_profile}")
                
                with col2:
                    st.markdown("**峰值功率 (入射)**")
                    st.metric("P₁ 峰值", f"{results['P1_peak_in']/1e3:.2f} kW")
                    if process_type == 'SFG':
                        st.metric("P₂ 峰值", f"{results['P2_peak_in']/1e3:.2f} kW")
                
                with col3:
                    st.markdown("**单脉冲能量**")
                    st.metric("E₁ 入", f"{results['E1_in']*1e6:.2f} µJ")
                    st.metric("E₃ 出", f"{results['E3_out']*1e6:.2f} µJ")
                    if results['E1_in'] > 0:
                        eff = results['E3_out'] / results['E1_in'] * 100
                        st.metric("转换效率", f"{eff:.2f}%")
                
                # 脉冲能量演化
                st.subheader("📈 脉冲能量演化")
                fig = sim.plot_pulse_energy_evolution(results)
                st.pyplot(fig)
                plt.close(fig)
                
    except Exception as e:
        st.error(f"❌ 仿真出错: {str(e)}")
        st.exception(e)

# ============ 连续光结果显示（只有运行成功后才显示） ============
if 'cw_results' in st.session_state and st.session_state.get('cw_mode', False) and st.session_state.get('run_success', False):
    import matplotlib.pyplot as plt
    
    results = st.session_state['cw_results']
    sim = st.session_state['sim']
    L_mm_cached = st.session_state['L_mm']
    process_type_cached = st.session_state['process_type']
    beam_profile_code_cached = st.session_state['beam_profile_code']
    lam3_nm_cached = st.session_state['lam3_nm']
    beam_profile_cached = st.session_state['beam_profile']
    P1_peak_cached = st.session_state['P1_peak']
    P2_peak_cached = st.session_state['P2_peak']
    P3_peak_cached = st.session_state['P3_peak']
    
    z = results['z_axis']
    P1 = results['P1_evol']
    P2 = results['P2_evol']
    P3 = results['P3_evol']
    P1_out = results['P1_out']
    P2_out = results['P2_out']
    P3_out = results['P3_out']
    
    # 对于高斯光斑，获取积分后的总功率演化
    if beam_profile_code_cached == 'gaussian' and 'P1_z' in results:
        P1_plot = results['P1_z']  # 总功率演化
        P2_plot = results['P2_z']
        P3_plot = results['P3_z']
        ylabel_suffix = "(Total Power)"
    else:
        P1_plot = P1  # 平顶光斑直接是功率
        P2_plot = P2
        P3_plot = P3
        ylabel_suffix = ""
    
    # 显示结果
    st.success(f"✅ {process_type_cached} 仿真完成！（{beam_profile_cached}）")
    
    # 计算功率密度和转换效率分析
    w0 = results['w0']
    w0_mm = w0 * 1e3
    if beam_profile_code_cached == 'gaussian':
        # 高斯光斑峰值功率密度
        I1_peak = results['I1_peak_in']
        I1_peak_MW_cm2 = I1_peak / 1e10  # W/m² -> MW/cm²
        conversion_eff = P3_out / P1_peak_cached * 100 if process_type_cached == 'SHG' else P3_out / (P1_peak_cached + P2_peak_cached) * 100
        
        # 判断转换状态
        if conversion_eff < 5:
            regime_msg = "⚠️ **低转换区**（<5%）：输出光斑形状与输入几乎相同"
            regime_color = "blue"
        elif conversion_eff < 20:
            regime_msg = "📊 **中等转换区**（5-20%）：开始出现轻微光斑畸变"
            regime_color = "orange"
        else:
            regime_msg = "🔥 **高转换/泵浦耗尽区**（>20%）：明显的光斑形状变化！"
            regime_color = "red"
        
        st.markdown(f":{regime_color}[{regime_msg}]")
        st.caption(f"峰值功率密度：{I1_peak_MW_cm2:.3f} MW/cm² | 光斑半径 w₀ = {w0_mm:.3f} mm")
    
    process_label = "倍频" if process_type_cached == "SHG" else "和频"
    
    if process_type_cached == 'SHG':
        # SHG：只显示P1和P3
        col1, col2 = st.columns(2)
        with col1:
            st.metric(f"{process_label}波长 λ₃", f"{lam3_nm_cached:.2f} nm")
            st.metric("P₁ 基频光出射", f"{P1_out:.4f} W", f"{(P1_out-P1_peak_cached)/P1_peak_cached*100:.2f}%" if P1_peak_cached > 0 else "N/A")
        
        with col2:
            st.metric("P₃ 倍频光出射", f"{P3_out:.4f} W")
            conversion_eff = P3_out / P1_peak_cached * 100 if P1_peak_cached > 0 else 0
            st.metric("转换效率", f"{conversion_eff:.2f}%")
            energy_total = P1_out + P3_out
            energy_in = P1_peak_cached
            st.metric("能量守恒检查", f"{energy_total:.4f} W", f"{(energy_total-energy_in)/energy_in*100:.2f}%" if energy_in > 0 else "N/A")
    else:
        # SFG：显示P1, P2, P3
        col1, col2 = st.columns(2)
        with col1:
            st.metric(f"{process_label}波长 λ₃", f"{lam3_nm_cached:.2f} nm")
            st.metric("P₁ 出射功率", f"{P1_out:.4f} W", f"{(P1_out-P1_peak_cached)/P1_peak_cached*100:.2f}%" if P1_peak_cached > 0 else "N/A")
            st.metric("P₂ 出射功率", f"{P2_out:.4f} W", f"{(P2_out-P2_peak_cached)/P2_peak_cached*100:.2f}%" if P2_peak_cached > 0 else "N/A")
        
        with col2:
            st.metric("P₃ 出射功率", f"{P3_out:.4f} W")
            conversion_eff = P3_out / (P1_peak_cached + P2_peak_cached) * 100 if (P1_peak_cached + P2_peak_cached) > 0 else 0
            st.metric("转换效率", f"{conversion_eff:.2f}%")
            energy_total = P1_out + P2_out + P3_out
            energy_in = P1_peak_cached + P2_peak_cached + P3_peak_cached
            st.metric("能量守恒检查", f"{energy_total:.4f} W", f"{(energy_total-energy_in)/energy_in*100:.2f}%" if energy_in > 0 else "N/A")
    
    # 绘制功率演化曲线
    st.subheader("📈 功率演化曲线")
    if beam_profile_code_cached == 'gaussian':
        st.caption("（显示空间积分后的总功率）")
    fig = sim.plot_results(z, P1_plot, P2_plot, P3_plot)
    st.pyplot(fig)
    plt.close(fig)
    
    # 光斑热力图 - 3行2列显示三束光的输入输出对比
    st.subheader("🔥 三束光光斑演化对比")
    st.write("左列：晶体入口 (z=0) | 右列：晶体出口 (z=L)")
    
    # 获取参数
    z_axis = results['z_axis']
    r_axis = results['r_axis']
    w0_local = results['w0']
    
    # 获取输入和输出的强度分布
    I1_in = results['I1_2d'][0, :]
    I1_out = results['I1_2d'][-1, :]
    I2_in = results['I2_2d'][0, :]
    I2_out = results['I2_2d'][-1, :]
    I3_in = results['I3_2d'][0, :]
    I3_out = results['I3_2d'][-1, :]
    
    # 创建2D圆形网格
    r_mm = r_axis * 1e3
    n_points = 150
    x = np.linspace(-r_mm[-1], r_mm[-1], n_points)
    y = np.linspace(-r_mm[-1], r_mm[-1], n_points)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    
    from scipy.interpolate import interp1d
    
    # 插值函数
    def radial_to_2d(I_r, r_mm, R):
        interp_func = interp1d(r_mm, I_r, kind='linear', bounds_error=False, fill_value=0)
        return interp_func(R) / 1e4  # 转换为 W/cm²
    
    # 转换为2D分布
    I1_in_2d = radial_to_2d(I1_in, r_mm, R)
    I1_out_2d = radial_to_2d(I1_out, r_mm, R)
    I2_in_2d = radial_to_2d(I2_in, r_mm, R)
    I2_out_2d = radial_to_2d(I2_out, r_mm, R)
    I3_in_2d = radial_to_2d(I3_in, r_mm, R)
    I3_out_2d = radial_to_2d(I3_out, r_mm, R)
    
    # 确定每束光的色标范围（输入输出统一）
    vmax1 = max(np.max(I1_in_2d), np.max(I1_out_2d)) if max(np.max(I1_in_2d), np.max(I1_out_2d)) > 0 else 1
    vmax2 = max(np.max(I2_in_2d), np.max(I2_out_2d)) if max(np.max(I2_in_2d), np.max(I2_out_2d)) > 0 else 1
    vmax3 = max(np.max(I3_in_2d), np.max(I3_out_2d)) if max(np.max(I3_in_2d), np.max(I3_out_2d)) > 0 else 1
    
    # 波长和标签
    lam1_nm = sim.lam1 * 1e9
    lam2_nm = sim.lam2 * 1e9
    lam3_nm = sim.lam3 * 1e9
    
    if process_type_cached == 'SHG':
        labels = [f'P1 Fundamental ({lam1_nm:.1f}nm)', f'P2 (=P1)', f'P3 SH ({lam3_nm:.1f}nm)']
    else:
        labels = [f'P1 Signal ({lam1_nm:.1f}nm)', f'P2 Idler ({lam2_nm:.1f}nm)', f'P3 Sum ({lam3_nm:.1f}nm)']
    
    # 创建3行2列的图
    fig_beams, axes = plt.subplots(3, 2, figsize=(12, 15))
    w0_mm_local = w0_local * 1e3
    
    data_pairs = [
        (I1_in_2d, I1_out_2d, vmax1, labels[0]),
        (I2_in_2d, I2_out_2d, vmax2, labels[1]),
        (I3_in_2d, I3_out_2d, vmax3, labels[2]),
    ]
    
    for row, (I_in, I_out, vmax, label) in enumerate(data_pairs):
        # 输入光斑
        ax_in = axes[row, 0]
        im_in = ax_in.pcolormesh(X, Y, I_in, cmap='jet', shading='auto', vmin=0, vmax=vmax)
        plt.colorbar(im_in, ax=ax_in, label='I (W/cm²)')
        circle_in = plt.Circle((0, 0), w0_mm_local, fill=False, color='white', linestyle='--', linewidth=1.5)
        ax_in.add_patch(circle_in)
        ax_in.set_xlabel('x (mm)')
        ax_in.set_ylabel('y (mm)')
        ax_in.set_title(f'{label}\nInput (z=0) | Peak: {np.max(I_in):.2e} W/cm²')
        ax_in.set_aspect('equal')
        
        # 输出光斑
        ax_out = axes[row, 1]
        im_out = ax_out.pcolormesh(X, Y, I_out, cmap='jet', shading='auto', vmin=0, vmax=vmax)
        plt.colorbar(im_out, ax=ax_out, label='I (W/cm²)')
        circle_out = plt.Circle((0, 0), w0_mm_local, fill=False, color='white', linestyle='--', linewidth=1.5)
        ax_out.add_patch(circle_out)
        ax_out.set_xlabel('x (mm)')
        ax_out.set_ylabel('y (mm)')
        ax_out.set_title(f'{label}\nOutput (z=L) | Peak: {np.max(I_out):.2e} W/cm²')
        ax_out.set_aspect('equal')
    
    fig_beams.suptitle(f'CW Beam Profile Evolution ({beam_profile_cached})', fontsize=14, y=1.01)
    fig_beams.tight_layout()
    st.pyplot(fig_beams)
    plt.close(fig_beams)
    
    # 显示峰值强度汇总信息
    st.markdown("---")
    st.markdown("**峰值强度汇总 (W/cm²)**")
    col_sum1, col_sum2, col_sum3 = st.columns(3)
    with col_sum1:
        st.write(f"**{labels[0]}**")
        st.write(f"入射: {np.max(I1_in_2d):.2e}")
        st.write(f"出射: {np.max(I1_out_2d):.2e}")
    with col_sum2:
        st.write(f"**{labels[1]}**")
        st.write(f"入射: {np.max(I2_in_2d):.2e}")
        st.write(f"出射: {np.max(I2_out_2d):.2e}")
    with col_sum3:
        st.write(f"**{labels[2]}**")
        st.write(f"入射: {np.max(I3_in_2d):.2e}")
        st.write(f"出射: {np.max(I3_out_2d):.2e}")
    
    # 数据表格
    with st.expander("📋 查看详细数据"):
        import pandas as pd
        if beam_profile_code_cached == 'gaussian':
            # 高斯光斑：同时显示中心强度和总功率
            df = pd.DataFrame({
                'z (mm)': z * 1e3,
                'P1 total (W)': P1_plot,
                'P2 total (W)': P2_plot,
                'P3 total (W)': P3_plot,
                'I1 center (W/m²)': P1,
                'I2 center (W/m²)': P2,
                'I3 center (W/m²)': P3
            })
        else:
            df = pd.DataFrame({
                'z (mm)': z * 1e3,
                'P1 (W)': P1,
                'P2 (W)': P2,
                'P3 (W)': P3
            })
        st.dataframe(df, height=300)

# ============ 脉冲空间分析模式的结果显示 ============
if 'pulse_spatial_results' in st.session_state and st.session_state.get('pulse_mode') == 'spatial' and st.session_state.get('run_success', False):
    import matplotlib.pyplot as plt
    import numpy as np
    
    results = st.session_state['pulse_spatial_results']
    sim = st.session_state['sim']
    process_type_cached = st.session_state.get('pulse_process_type', 'SFG')
    beam_profile_cached = st.session_state['pulse_beam_profile']
    
    st.markdown("---")
    st.subheader("🌡️ 三束光脉冲能量密度演化对比 (Fluence)")
    st.write("左列：晶体入口 (z=0) | 右列：晶体出口 (z=L)")
    st.write("*将所有时间切片的光斑叠加后得到的脉冲能量密度分布*")
    
    # 获取参数
    r_axis = results['r_axis']
    w0_local = results['w0']
    
    # 获取输入和输出的能量密度分布
    F1_in = results['F1_in']
    F1_out = results['F1_out']
    F2_in = results['F2_in']
    F2_out = results['F2_out']
    F3_in = results['F3_in']
    F3_out = results['F3_out']
    
    # 创建2D圆形网格
    r_mm = r_axis * 1e3
    n_points = 150
    x = np.linspace(-r_mm[-1], r_mm[-1], n_points)
    y = np.linspace(-r_mm[-1], r_mm[-1], n_points)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    
    from scipy.interpolate import interp1d
    
    # 插值函数
    def radial_to_2d(F_r, r_mm, R):
        interp_func = interp1d(r_mm, F_r, kind='linear', bounds_error=False, fill_value=0)
        return interp_func(R) / 1e4  # 转换为 J/cm²
    
    # 转换为2D分布
    F1_in_2d = radial_to_2d(F1_in, r_mm, R)
    F1_out_2d = radial_to_2d(F1_out, r_mm, R)
    F2_in_2d = radial_to_2d(F2_in, r_mm, R)
    F2_out_2d = radial_to_2d(F2_out, r_mm, R)
    F3_in_2d = radial_to_2d(F3_in, r_mm, R)
    F3_out_2d = radial_to_2d(F3_out, r_mm, R)
    
    # 确定每束光的色标范围（输入输出统一）
    vmax1 = max(np.max(F1_in_2d), np.max(F1_out_2d)) if max(np.max(F1_in_2d), np.max(F1_out_2d)) > 0 else 1e-10
    vmax2 = max(np.max(F2_in_2d), np.max(F2_out_2d)) if max(np.max(F2_in_2d), np.max(F2_out_2d)) > 0 else 1e-10
    vmax3 = max(np.max(F3_in_2d), np.max(F3_out_2d)) if max(np.max(F3_in_2d), np.max(F3_out_2d)) > 0 else 1e-10
    
    # 波长和标签（使用缓存的数据）
    lam1_nm = st.session_state.get('pulse_lam1_nm', 1064.0)
    lam2_nm = st.session_state.get('pulse_lam2_nm', 532.0)
    lam3_nm = st.session_state.get('pulse_lam3_nm', 355.0)
    
    if process_type_cached == 'SHG':
        labels = [f'P1 Fundamental ({lam1_nm:.1f}nm)', f'P2 (=P1)', f'P3 SH ({lam3_nm:.1f}nm)']
    else:
        labels = [f'P1 Signal ({lam1_nm:.1f}nm)', f'P2 Idler ({lam2_nm:.1f}nm)', f'P3 Sum ({lam3_nm:.1f}nm)']
    
    # 创建3行2列的图
    fig_beams, axes = plt.subplots(3, 2, figsize=(12, 15))
    w0_mm_local = w0_local * 1e3
    
    data_pairs = [
        (F1_in_2d, F1_out_2d, vmax1, labels[0]),
        (F2_in_2d, F2_out_2d, vmax2, labels[1]),
        (F3_in_2d, F3_out_2d, vmax3, labels[2]),
    ]
    
    for row, (F_in, F_out, vmax, label) in enumerate(data_pairs):
        # 输入光斑
        ax_in = axes[row, 0]
        im_in = ax_in.pcolormesh(X, Y, F_in, cmap='jet', shading='auto', vmin=0, vmax=vmax)
        plt.colorbar(im_in, ax=ax_in, label='F (J/cm²)')
        circle_in = plt.Circle((0, 0), w0_mm_local, fill=False, color='white', linestyle='--', linewidth=1.5)
        ax_in.add_patch(circle_in)
        ax_in.set_xlabel('x (mm)')
        ax_in.set_ylabel('y (mm)')
        ax_in.set_title(f'{label}\nInput (z=0) | Peak: {np.max(F_in):.2e} J/cm²')
        ax_in.set_aspect('equal')
        
        # 输出光斑
        ax_out = axes[row, 1]
        im_out = ax_out.pcolormesh(X, Y, F_out, cmap='jet', shading='auto', vmin=0, vmax=vmax)
        plt.colorbar(im_out, ax=ax_out, label='F (J/cm²)')
        circle_out = plt.Circle((0, 0), w0_mm_local, fill=False, color='white', linestyle='--', linewidth=1.5)
        ax_out.add_patch(circle_out)
        ax_out.set_xlabel('x (mm)')
        ax_out.set_ylabel('y (mm)')
        ax_out.set_title(f'{label}\nOutput (z=L) | Peak: {np.max(F_out):.2e} J/cm²')
        ax_out.set_aspect('equal')
    
    fig_beams.suptitle(f'Pulsed Beam Fluence Evolution ({beam_profile_cached})', fontsize=14, y=1.01)
    fig_beams.tight_layout()
    st.pyplot(fig_beams)
    plt.close(fig_beams)
    
    # 显示峰值能量密度汇总信息
    st.markdown("---")
    st.markdown("**峰值能量密度汇总 (J/cm²)**")
    col_sum1, col_sum2, col_sum3 = st.columns(3)
    with col_sum1:
        st.write(f"**{labels[0]}**")
        st.write(f"入射: {np.max(F1_in_2d):.2e}")
        st.write(f"出射: {np.max(F1_out_2d):.2e}")
    with col_sum2:
        st.write(f"**{labels[1]}**")
        st.write(f"入射: {np.max(F2_in_2d):.2e}")
        st.write(f"出射: {np.max(F2_out_2d):.2e}")
    with col_sum3:
        st.write(f"**{labels[2]}**")
        st.write(f"入射: {np.max(F3_in_2d):.2e}")
        st.write(f"出射: {np.max(F3_out_2d):.2e}")
    
    # 能量密度数据
    with st.expander("📋 能量密度径向分布数据"):
        import pandas as pd
        
        if process_type_cached == 'SHG':
            df_fluence = pd.DataFrame({
                'r (mm)': r_mm,
                'F1_in (J/cm²)': F1_in / 1e4,
                'F3_in (J/cm²)': F3_in / 1e4,
                'F1_out (J/cm²)': F1_out / 1e4,
                'F3_out (J/cm²)': F3_out / 1e4,
            })
        else:
            df_fluence = pd.DataFrame({
                'r (mm)': r_mm,
                'F1_in (J/cm²)': F1_in / 1e4,
                'F2_in (J/cm²)': F2_in / 1e4,
                'F3_in (J/cm²)': F3_in / 1e4,
                'F1_out (J/cm²)': F1_out / 1e4,
                'F2_out (J/cm²)': F2_out / 1e4,
                'F3_out (J/cm²)': F3_out / 1e4,
            })
        st.dataframe(df_fluence, height=300)

# ============ 默认页面显示（当没有运行结果时） ============
if not ('cw_results' in st.session_state or 'pulse_spatial_results' in st.session_state):
    st.info("👈 请在左侧设置参数，然后点击\"运行仿真\"按钮")
    
    # 显示理论说明
    st.subheader("📖 理论背景")
    st.markdown("""
    ### 耦合波方程（Boyd体系）
    
    本程序支持两种非线性过程：
    
    #### 1. 和频产生 (SFG): ω₃ = ω₁ + ω₂
    
    耦合波方程：
    
    $$\\frac{dA_1}{dz} = iK_1 A_3 A_2^* e^{-i\\Delta k z} - \\frac{\\alpha_1}{2}A_1$$
    
    $$\\frac{dA_2}{dz} = iK_2 A_3 A_1^* e^{-i\\Delta k z} - \\frac{\\alpha_2}{2}A_2$$
    
    $$\\frac{dA_3}{dz} = iK_3 A_1 A_2 e^{i\\Delta k z} - \\frac{\\alpha_3}{2}A_3$$
    
    其中耦合系数：$K_i = \\frac{2\\omega_i d_{eff}}{n_i c}$ （**有系数2**）
    
    #### 2. 二次谐波产生 (SHG): ω₃ = 2ω₁
    
    SHG是SFG的特殊情况（ω₁ = ω₂），但耦合系数**没有系数2**：
    
    $$K_i = \\frac{\\omega_i d_{eff}}{n_i c}$$ （**无系数2**）
    
    #### 脉冲光仿真
    
    对于纳秒脉冲激光，采用**准静态近似**（切片法）：
    - 将脉冲切成多个时间片
    - 每个时间片独立求解耦合波方程
    - 重组得到输出脉冲形状
    
    ### 参数说明
    - 光强定义：$I = 2n\\epsilon_0 c|A|^2$ (Boyd体系)
    - $\\Delta k = k_3 - k_1 - k_2$ 为相位失配
    - $\\alpha_i$ 为吸收系数
    
    ---
    
    ### ⚡ 转换效率与光斑形状变化
    
    #### 💡 为什么转换效率与功率密度平方成正比？
    
    **小信号近似（低转换效率 <5%）**：
    
    对于SHG，转换效率为：
    
    $$\\eta \\propto I_1^2 L^2 \\propto \\left(\\frac{P}{\\pi w_0^2}\\right)^2 L^2$$
    
    关键结论：
    - **功率翻倍** → 转换效率提高 **4倍**
    - **光斑半径减半** → 功率密度提高 4倍 → 转换效率提高 **16倍**！
    
    #### 🔥 光斑形状什么时候会变化？
    
    | 转换效率 | 转换状态 | 光斑形状变化 |
    |---------|---------|-------------|
    | < 5% | 小信号区 | ❌ 几乎无变化（各点独立转换） |
    | 5-20% | 中等转换 | ⚠️ 开始出现轻微畸变 |
    | > 20% | 泵浦耗尽 | ✅ **明显变化**：泵浦变平，倍频变尖 |
    
    **物理原因**：
    - 光斑中心强度最高 → 转换最多 → 泵浦消耗最多
    - 边缘强度低 → 转换少 → 相对保留更多
    - 结果：P1输出变平，P3输出更尖锐
    
    """)

