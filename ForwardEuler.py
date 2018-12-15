'''
Forward method to solve 1D reaction-diffusion equation:
        u_t = k * u_xx
    
with Neumann boundary conditions 
at x=0: u_x(0,t) = 0 = sin(2*np.pi)
at x=L: u_x(L,t) = 0 = sin(2*np.pi)

with L = 1 and initial conditions:
u(x,0) = (1.0/2.0)+ np.cos(2.0*np.pi*x) - (1.0/2.0)*np.cos(3*np.pi*x)

u_x(x,t) = (-4.0*(np.pi**2))np.exp(-4.0*(np.pi**2)*t)*np.cos(2.0*np.pi*x) + 
            (9.0/2.0)*(np.pi**2)*np.exp(-9.0*(np.pi**2)*t)*np.cos(3*np.pi*x))
'''

def ForwardEuler(M, lambd, T = 0.5, L = 1, k = 1): 
    #Parameters needed to solve the equation within the explicit method
    #M = GRID POINTS on space interval
    N = (M**2) #GRID POINTS on time interval

    # ---- Length of the wire in x direction ---- 
    x0 = 0
    xL = L

    # ----- Spatial discretization step -----
    dx = (xL - x0)/(M-1)

    # ---- Final time ---- 
    t0 = 0
    tF = T 

    # ----- Time step -----
    dt = (tF - t0)/(N-1)

    #lambd = dt*k/dx**2
    
    # ----- Creates grids -----
    xspan = np.linspace(x0, xL, M)
    tspan = np.linspace(t0, tF, N)

    # ----- Initializes matrix solution U -----
    U = np.zeros((M, N))

    # ----- Initial condition -----
    U[:,0] = (1.0/2.0)+ np.cos(2.0*np.pi*xspan) - (1.0/2.0)*np.cos(3*np.pi*xspan)

    # ----- Neumann boundary conditions -----
    """
    To implement these boundary conditions, we again use “false points”, x_0 and x_N+1 which are external points. 
    We use a difference to approximate ∂u/∂x (xL,t) and set it equal to the desired boundary condition:
    
    BEFORE: 
    U[0,:] = 0.0
    U[-1,:] = 0.0
    
    NOW:
    U[0,:] = (-3*U[0,:] + 4*U[1,:] - U[2,:])/2*dx = 0
    U[-1,:] = (-3*U[-1,:] + 4*U[-2,:] - U[-3,:])/2*dx = 0
    """
    
    f = np.arange(1, N+1)
    f = (-3*U[0,:] + 4*U[1,:] - U[2,:])/2*dx
    U[0,:] = (4*U[1,:] - U[2,:])/3
    
    g = np.arange(1, N+1)
    g = (-3*U[-1,:] + 4*U[-2,:] - U[-3,:])/2*dx
    U[-1,:] = (4*U[-2,:] - U[-3,:])/3
    
    # -----  ftcs -----
    for k in range(0, N-1):
        for i in range(1, M-1):
            U[i, k+1] = lambd*U[i-1, k] + (1-2*lambd)*U[i,k] + lambd*U[i+1,k] 
    
    return (U, tspan, xspan)

U, tspan, xspan = ForwardEuler(M = 14, lambd = 1.0/6.0)
Uexact, x, t = ExactSolution(M = 14)
surfaceplot(U, Uexact, tspan, xspan, M = 14)
