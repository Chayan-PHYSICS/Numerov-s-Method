import  numpy as np
import matplotlib.pyplot as plt
from scipy import special

n = int(input("How many state you want  = "))

def integration(y,x):
    '''y: Integrand : x dependent
       x: Independent variable on which y depends
  return: Integration of y
    '''
    return(np.abs(x[1] - x[0])/3) * (y[0] + 4* np.sum(y[1: len(y)-1: 2]) + 2 * np.sum(y[2:len(y)-2 : 2]) + y[-1])

total_psi_dict = {}
for j in range(n+1):
    E_max = 20.0
    E_min = 0
    dx = 0.001# spatile grid point
    x_grid_max = 11
    energy_guess = (E_max + E_min)/2

    #special grid point 
    x_arr = np.arange(0,x_grid_max,dx)
    PSI = np.zeros(len(x_arr)) #psi array for wavefunction

    # potential energy function
    V = 0.5*1*x_arr**2

    def Classical_inversion_point(): 
        """Return ==> (1) Classical inversion point\n
        (2)Array from zero to classical inversion point \n
        (3)Array from x_grid_max to classical inversion point"""
        #set energy
        energy_array = energy_guess*np.ones(len(x_arr))
        # find classical inversion point
        I_point = 0
        energy_diff = energy_array - V
        
        for i in range (len(energy_diff)-1):
            prd = energy_diff[i]*energy_diff[i-1]
            if energy_diff[i] == 0.0:
                energy_diff[i] = 10**(-10)
            if prd < 0 :
                I_point = i
            
        # Array from 0 to classical inversion point
        arr1_cls = x_arr[:I_point+1]
        # Array from clasical invertion point to end point
        arr2_cls = x_arr[I_point:][::-1]

        return arr1_cls,arr2_cls,I_point

    arr1, arr2, I_point  = Classical_inversion_point()

    def Wave_function():
        """ Main function to creat Wavefuncton of a given energy\n
            Return ==> (1) Psi (Wavefunction)\n
                    (2) Discontinuty in first order derivative of wavefunction"""
        y = 1 - ((dx)**2/12)*(x_arr**2-2*energy_guess)
        y_in = y[I_point:][::-1]
        # outward integration
        # for odd wavefunction
        if j/2 != int(j/2):  
            PSI[0] = 0 
            PSI[1] = dx
        # For even function
        else:
            PSI[0] = 1
            PSI[1] = 0.5*(12 - 10*y[0])*PSI[0]/y[1]
        # Numerove methode
        for i in range(1, len(arr1)-1):
            PSI[i+1] = ((12 - 10*y[i])*PSI[i]-PSI[i-1]*y[i-1])/y[i+1] 

        #Inward integration
        psi_in = np.zeros(len(x_arr[I_point:]))
        psi_in[0] = 0
        psi_in[1] = dx
        for i in range(1, len(arr2)-1):
            psi_in[i+1] = ((12 - 10*y_in[i])*psi_in[i]-psi_in[i-1]*y_in[i-1])/y_in[i+1]

        rescaling_factor = PSI[I_point]/psi_in[-1]
        psi_in_actuall = (psi_in*rescaling_factor )[::-1]

        for i in range(len(arr2)):
            PSI[I_point +i] = psi_in_actuall[i]

        #Discontinuty in the first order derivative0

        discontinuty_1st_Order = (PSI[I_point - 1] + PSI[I_point + 1]- (14 -12*y[I_point])*PSI[I_point])/dx
        return PSI, discontinuty_1st_Order

    while True:
        # arr1, arr2, I_point  = Classical_inversion_point()
        psi, discontinuty_1st_Order = Wave_function()
        
        
        #counting the number of times wave function change its sing
        count_n = 0
        for i in range (len(psi)-1):
            if psi[i] * psi[i+1] < 0:
                count_n = count_n + 1

        # modify the count_n for odd and even function
        if j/2 != int(j/2):
            count_n = 2*count_n + 1
        else:
            count_n = 2* count_n

        # modify energy according to count_n
        if count_n != j:
            if count_n > j :
                E_max = energy_guess
            elif count_n < j :
                E_min = energy_guess
            energy_guess = 0.5*(E_max + E_min)
            continue
        # Derivative matching
        if discontinuty_1st_Order* psi[I_point]>= 0:
            E_max = energy_guess
        else:
            E_min = energy_guess
        energy_guess = 0.5*(E_max + E_min)
        # Loop exiting condition
        if (E_max - E_min) <= 10 **(-10):
            break
    print(f'Energy{j} =', energy_guess)

    # Total psi (wavefunction)
    total_psi = []
    for i in PSI[-1:-len(PSI):-1]:
        total_psi.append((-1)**j*i)
    for i in PSI:
        total_psi.append(i)

    total_psi = np.array(total_psi)
    # Normalized total_psi(Wavefunction)
    total_psi = total_psi * np.sqrt(1/ integration((total_psi)**2, np.linspace(-x_grid_max,+x_grid_max,len(total_psi)))) 

    # Add key value in the Psi dictionary
    key = 'total_psi'+ str(j)
    total_psi_dict[key] = total_psi

# Checking the orthogonality of the wavefunction
if input("Are you want to check orthonormality?(y/n): ") == 'y':
    print()
    n1,n2 = input(f'Enter two value between 0 & {n}: ').split(',') #with comma
    total_xarray_1 = np.linspace(-x_grid_max,+x_grid_max,len(total_psi_dict['total_psi'+ str(n1)]))
    if n1 == n2 :
        orthonormality = integration((total_psi_dict['total_psi'+ str(n1)]) ** 2, total_xarray_1)
        print("the wave function is normalized with value = ", orthonormality)
    else:
        orthonormality = integration((total_psi_dict['total_psi'+ str(n1)]) * total_psi_dict['total_psi'+ str(n2)], total_xarray_1)
        print("The wavefunction are orthogonal with value = ",orthonormality)
else:
    pass
print()
# Expectation values & Uncertainty 
if input("Are you want to check position & momentum uncertainty?(y/n): ") == 'y':
    #Expectation value of position on a wavefunction
    state_number = input(f'Enter a value between 0 & {n}: ')
    total_xarray_2 = np.linspace(-x_grid_max,+x_grid_max,len(total_psi_dict['total_psi'+ str(state_number)]))
    expt_x  = integration(total_psi_dict['total_psi'+ str(state_number)]**2 * total_xarray_2 ,total_xarray_2) 
    expt_x2 = integration(total_psi_dict['total_psi'+ str(state_number)]**2 * total_xarray_2**2 ,total_xarray_2)
    del_x = (expt_x2 - expt_x**2)**0.5
    # print("expt_x = ",expt_x)

    # Expectation value of momentum on a wavefunction
    dpsi_dx = np.gradient(total_psi_dict['total_psi'+ str(state_number)], total_xarray_2)
    d2psi_dx2 = np.gradient(dpsi_dx, total_xarray_2)
    expt_p = integration(total_psi_dict['total_psi'+ str(state_number)] * dpsi_dx, total_xarray_2)
    expt_p2 = -integration(total_psi_dict['total_psi'+ str(state_number)] * d2psi_dx2, total_xarray_2)
    del_p = (expt_p2 + expt_p **2) ** 0.5

    print(del_x)
    print(del_p)
    print(del_p * del_x)
    print(del_x * del_p >= 0.49999) 
else:
    pass

"""Ploting the wavefunction"""
print()
if input("Are you want to see the wavefunction?(y/n): ") == 'y':
    
    # plt.style.use('ggplot')
    plt.rc('font', **{'family': 'serif', 'size':10})
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax2 = fig.add_subplot(122)

    wave_N = str(input("what wavefunction you want to print: " ))

    total_xarray = np.linspace(-x_grid_max,+x_grid_max,len(total_psi_dict['total_psi'+ str(wave_N)]))

    def actual_wavefunction(wave_N):
        hermite = special.hermite(wave_N)(total_xarray)
        return  (1/np.pi **0.25) * 1 /(2**wave_N * np.math.factorial(wave_N))**0.5 *\
                hermite*np.exp(-total_xarray **2/ 2)

    if int(wave_N )== 0:
        plt.suptitle(r'GROUND STATE')
    else:
        plt.suptitle(f'{wave_N} EXCITED STATE')

    ax1.plot(total_xarray,actual_wavefunction(int(wave_N)),'b',label = 'Analitycal Solution')
    ax1.plot(total_xarray,total_psi_dict['total_psi'+ str(wave_N)], '--r', label='Numarical Solution ')
    # first_state = np.gradient(total_psi_dict['total_psi'+ str(wave_N)], total_xarray)
    # ax2.plot(first_state,total_xarray)
    

    print()
    print("from analytic solution")
    expt_x_2  = integration(actual_wavefunction(int(wave_N))**2 * total_xarray ,total_xarray)
    expt_x2_2 = integration(actual_wavefunction(int(wave_N))**2 * total_xarray**2 ,total_xarray)
    del_x_2 = (expt_x2_2 - expt_x_2**2)**0.5
    print("del_x =",del_x_2)
 
    dpsi_dx_2= np.gradient(actual_wavefunction(int(wave_N)), total_xarray)
    d2psi_dx2_2 = np.gradient(dpsi_dx_2, total_xarray)
    expt_p_2 = integration(actual_wavefunction(int(wave_N)) * dpsi_dx_2, total_xarray)
    expt_p2_2 = -integration(actual_wavefunction(int(wave_N)) * d2psi_dx2_2, total_xarray)
    del_p_2 = (expt_p2_2 + expt_p_2 **2) ** 0.5
    print("del_p=",del_p_2)

    # ax1.set_xticks(np.linspace(-x_grid_max , +x_grid_max , 10))
    # ax2.set_xticks(np.linspace(-x_grid_max , +x_grid_max , 10))
    # ax1.grid(True, linewidth=0.7, alpha=0.8)
    # ax1.tick_params(color='red', width=5)
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.spines['left'].set_visible(False)
    # ax1.spines['bottom'].set_visible(False)
    plt.legend()
    plt.show ()
else:
    pass