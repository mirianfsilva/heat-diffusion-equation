# ---- Surface plot ----
def surfaceplot(U, Uexact, tspan, xspan, M): 
    N = M**2
    #meshgrid : Return coordinate matrices from coordinate vectors
    X, T = np.meshgrid(tspan, xspan)
    fig = plt.figure(figsize=plt.figaspect(0.5))
    #fig2 = plt.figure(figsize=plt.figaspect(0.5))
    #fig3 = plt.figure(figsize=plt.figaspect(0.5))
    
    # ---- Method Solution
    ax1 = fig.add_subplot(1, 3, 1,projection='3d')
    surf = ax1.plot_surface(X, T, U, linewidth=0, cmap=cm.jet, antialiased=True)
    ax1.set_title('Approximation')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Space')
    ax1.set_zlabel('Temp')
    
    
    # ---- Absolute Error ----
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    surf2 = ax2.plot_surface(X, T, abs(Uexact-U), linewidth=0, cmap=cm.jet, antialiased=True)
    ax2.set_title('Absolute Error')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Space')
    ax2.set_zlabel('Error')
    
    # ---- Relative Error ---- 
    Urelative = np.zeros((M, N))

    for i in range(0, M):
        for j in range(0, N):
              Urelative[i][j] = Uexact[i][j] - U[i][j] / Uexact[i][j]
    
          
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    surf3 = ax3.plot_surface(X, T, abs(Urelative), linewidth=0, cmap=cm.jet, antialiased=True)
    ax3.set_title('Relative Error')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Space')
    ax3.set_zlabel('Error')
    
    plt.tight_layout()
    ax1.view_init(30,230)
    ax2.view_init(30,230)
    ax3.view_init(30,230)
    
    plt.draw()
