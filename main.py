import matplotlib.pyplot as plt


def phi_next(dt, phi_i, u_i):
    return u_i * dt + phi_i


def x_next(dt, x_i, v_i):
    return v_i * dt + x_i


c = 0


def v_next(dt, x_i, phi_i, v_i, wq_x, k_x):
    return v_i - (wq_x * x_i + k_x * x_i * phi_i + c * v_i) * dt


def u_next(dt, x_i, phi_i, u_i, wq_phi, k_phi):
    return u_i - (wq_phi * phi_i + k_phi * x_i * phi_i + c * u_i) * dt


def calculate():
    with open("input.txt") as file:
        n, dt, tau, x0, phi0, v0, u0, k, m, I, D, a, b = map(float, file.readlines()[1].split())
    n = int(n)
    wq_phi = D / I
    k_phi = a / I
    wq_x = k / m
    k_x = b / m

    # Formalnosti
    t = 0
    T = []
    X = []
    PHI = []
    Ephi = []
    Ex = []
    E = []
    phi_i = phi0
    u_i = u0
    x_i = x0
    v_i = v0

    # Calculation
    while t <= tau:
        T.append(t)
        x = x_next(dt, x_i, v_i)
        v = v_next(dt=dt, x_i=x_i, phi_i=phi_i, v_i=v_i, wq_x=wq_x, k_x=k_x)
        phi = phi_next(dt, phi_i, u_i)
        u = u_next(dt=dt, x_i=x_i, phi_i=phi_i, u_i=u_i, wq_phi=wq_phi, k_phi=k_phi)
        X.append(round(x, n))
        Ex.append(round((k * x ** 2 / 2 + m * v ** 2 / 2), n))
        PHI.append(round(phi, n))
        Ephi.append(round(I * u ** 2 / 2 + D * phi**2 / 2, n))
        E.append(round(m * v ** 2 / 2 + k * x**2 / 2 + I * u ** 2 / 2 + D * phi**2 / 2, n))
        # Prepare to next step
        phi_i = phi
        u_i = u
        x_i = x
        v_i = v
        t += dt

    result = [T, X, PHI, Ex, Ephi, E]
    return result


def draw(t_x_phi):
    t = t_x_phi[0]
    x = t_x_phi[1]
    phi = t_x_phi[2]
    plt.plot(t, phi, 'b-')
    plt.plot(t, x, 'r-')
    plt.xlabel('t, c')
    plt.ylabel('x, phi')
    plt.grid()
    plt.show()


def draw1(t_x_phi):
    t = t_x_phi[0]
    e_x = t_x_phi[3]
    e_phi = t_x_phi[4]
    e = t_x_phi[5]
    plt.plot(t, e_x, 'b-')
    plt.plot(t, e_phi, 'r-')
    plt.plot(t, e, 'g-')
    plt.grid()
    plt.show()


def main():
    t_x_phi_array = calculate()
    draw(t_x_phi_array)
    draw1(t_x_phi_array)


if __name__ == "__main__":
    main()