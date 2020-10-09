import math
from numpy import array

def matrix_plot(matrix, a, b, mod_a, mod_b, n_a, n_b, density, N_tot):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
 
    print(matrix)
    print("\na and b vectors: ", a, b)
    print("Total length of a, b: %.2f %.2f"%(mod_a, mod_b))
    print("Density: %.2f pts/A^2.  Required: %2.f points"%(density, N_tot))
    print("Grid generated: ", n_a,"x",n_b,".  Total: %d points\n"%(n_a*n_b))
    
    plt.title("Projection on xy plane")
    plt.plot(matrix[:,0], matrix[:,1], 'o')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(matrix[:,0], matrix[:,1], matrix[:,2], c='r', marker='o', label="3D grid")
    x = [0,a[0],a[0]+b[0],b[0],0]
    y = [0,a[1],a[1]+b[1],b[1],0]
    z = [0,a[2],a[2]+b[2],b[2],0]
    ax.plot(x, y, z)
    plt.show()
    

def pes_grid(alats, density, mod=0, default=5):
    """
    alats must be a 2x3 matrix with the x, y, z components of a and b
    density is the number of points per unit A^2 required
    mod=0 construct a uniform grid respecting the ratio b/a
    mod=1 put default points along a and (b/a)*default along b
    default is the default value of points to take for a when mod=1
    """
    from numpy import zeros
    a = zeros(3) ; a = alats[0, :]
    b = zeros(3) ; b = alats[1, :]
    a_mod = math.sqrt(a[0]**2. + a[1]**2. + a[2]**2.)
    b_mod = math.sqrt(b[0]**2. + b[1]**2. + b[2]**2.)
    N_tot = round(density*a_mod*b_mod)
    
    ratio = b_mod/a_mod
    if mod == 0:
        n_a = int(round(math.sqrt(N_tot/ratio)))
        n_b = int(round(ratio*n_a))
    elif mod == 1:
        n_a = default
        n_b = int(round(ratio*default))
    k = 0
    dist_a_x = a[0]/n_a ; dist_a_y = a[1]/n_a ; dist_a_z = a[2]/n_a
    dist_b_x = b[0]/n_b ; dist_b_y = b[1]/n_b ; dist_b_z = b[2]/n_b
    matrix = zeros((n_a*n_b, 3))
    for i in range(0, n_a):
        for j in range(0, n_b):
            matrix[k, 0] = i*dist_a_x + j*dist_b_x
            matrix[k, 1] = i*dist_a_y + j*dist_b_y
            matrix[k, 2] = i*dist_a_z + j*dist_b_z
            k += 1
    matrix_plot(matrix, a, b, a_mod, b_mod, n_a, n_b, density, N_tot)

    return matrix
    
    
vectors = array([[3, 0, 0], 
                 [0.8, 4, 0.3]])

pes_grid(vectors, 3)
