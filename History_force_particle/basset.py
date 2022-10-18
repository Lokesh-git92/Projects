import numpy as np
from matplotlib import pyplot as plt

def particle_vel():
    #function that returns particle velcoity for a given time
    #uses drag+gravity in Maxey & Riley equation
    #implements forward difference finite difference method
    up = np.zeros(len(time))
    for it in range(1,len(time)):
        up[it] = (del_t/m_p) * (-6*np.pi*mu*0.5*dp*(uf[it-1]-up[it-1]) +((m_p-m_f)*g))

    return up

def particle_vel_history():
    #function that returns particle velcoity for a given time
    #uses drag+gravity in Maxey & Riley equation
    #impleemnts forward difference finite difference method
    up = np.zeros(len(time))
    g_t = np.zeros(len(time))
    C_b = 6 * (0.5*dp**2) * rho_f *np.sqrt(np.pi * nu)

    for it in range(1,len(time)):
        N=it-1
        start_ind = N- (N_win -1)
        #array to be susind in tail aprt computation
        Fi_old = np.zeros(10)
        Fi_new = np.zeros(10)
        if it == 1:
            F_h = (4/3)*C_b * g_t[N]*del_t

        elif N_win <= it:
            #for loop for calculating third term
            temp_sum = 0
            rev_ind = N-1
            for n in range(1,N-1):
                temp_sum += g_t[rev_ind] *(  ((n+(4/3)) / ((n+1)**1.5 +(n+1.5)*(n**0.5)) )\
                    + ((n-(4/3)) / ((n-1)**1.5 +(n-1.5)*(n**0.5) ))   )
                rev_ind -=1
            #history force term
            F_h = (4/3)*C_b * g_t[N]*del_t + C_b*g_t[0] * \
                ( ( N-(4/3) )/( (N-1)**1.5 +(N-1.5)* (N**0.5) )) \
                    +C_b * del_t * temp_sum

        else:
            #######################window part###########################
            #starting index
            start_ind = N- (N_win -1)
            #for loop for calculating third term
            temp_sum = 0
            #rev_index for calculating acceleration g_t
            rev_ind = N-1
            #third term summation
            for n in range(1,N-1):
                temp_sum += g_t[rev_ind] *(  ((n+(4/3)) / ((n+1)**1.5 +(n+1.5)*(n**0.5)) )\
                    + ((n-(4/3)) / ((n-1)**1.5 +(n-1.5)*(n**0.5) ))   )
                rev_ind -=1
            #history force
            F_h_win = (4/3)*C_b * g_t[N]*del_t + C_b*g_t[start_ind] * \
                ( ( N-(4/3) )/( (N-1)**1.5 +(N-1.5)* (N**0.5) )) \
                    +C_b * del_t * temp_sum
            F_h_tail = 0
            #############tail part###################
            for m in range(10):
                phi = lambda z : 1 + 0.5 *z + (1/6) * (z**2)
                #g_N
                first = g_t[start_ind] * (1- phi(-0.5*del_t/ti_bar[m]))
                #g_N+1
                second = g_t[start_ind - 1] * np.exp(-0.5*del_t/ti_bar[m]) * \
                    (phi(0.5*del_t/ti_bar[m]) -1)
                #t_win = N_w *del_t
                F_h_di = 2* C_b * (np.exp(1) *ti_bar[m]) * np.exp(-0.5 * N_win *del_t/ ti_bar[m]) \
                    * (first + second)
                F_h_re = np.exp(-0.5 * del_t/ti_bar[m]) * Fi_old[m]
                Fi_new[m] = F_h_di + F_h_re
                Fi_old[m] = Fi_new[m]
                
                F_h_tail += ai[m] * Fi_new[m]
            F_h = F_h_win + F_h_tail


        #particle velocity
        up[it] = (del_t/m_p) * (-6*np.pi*mu*0.5*dp*(uf[it-1]-up[it-1]) +((m_p-m_f)*g ) +F_h)
        #acceleration
        g_t[it] = ((uf[it] - uf[it-1])/del_t) - ((up[it] - up[it-1])/del_t)

    return up


if __name__ =="__main__":
    #print("Booyeah")
    #diameter, reynold numb of particle
    dp = 1.0
    Re_p = 0.1
    #densities
    rho_p = 1.1
    rho_f = 1.0
    #gravity
    g = 1.0
    #fluid viscosities
    mu = np.sqrt(rho_f*rho_p*(1-rho_f/rho_p)*g*(dp**3) /(18*Re_p) )
    nu = mu/ rho_f
    #particle's volume
    Vp=(np.pi * dp**3) / 6.0
    #mass of particle
    m_p = rho_p*Vp 
    #displaced mass of fluid by the presence of particle.
    m_f =  rho_f*Vp
    #particle relaxation time.
    tau_p=rho_p*dp**2/(18*mu) 
    N_tot=256
    time=np.linspace(0,0.1*np.pi,N_tot)
    del_t = time[1]-time[0]
    #quiescent fluid
    uf =np.zeros(len(time))
    #coeeficients for tail calculation
    ai=[0.23477481312586,0.28549576238194,0.28479416718255,0.26149775537574 \
        ,0.32056200511938,0.35354490689146,0.39635904496921,0.42253908596514 \
            ,0.48317384225265,0.63661146557001]
    ti_bar=[0.1,0.3,1,3,10,40,190,1000,6500,50000]
    #particle velocity calculated
    vel = particle_vel()
    #print(vel[-5:-1])
    N_win = 5 #python index from 0 to N-1 5points here from 0 to N-1
    vel_history = particle_vel_history()
    #print(vel_history[-5:-1])


    

    
    
    



    
     
