import matplotlib.pyplot as plt

from vector2 import vec2

grav_const = 6.67408e-11 # m3 kg-1 s-2

class body:
    def __init__(self, mass, radius):
        self.pos = vec2(0, 0)
        self.mass = mass
        self.radius = radius

class orbiter:
    def __init__(self, pos, vel, dropoff=2):
        self.pos = pos
        self.vel = vel
        self.dropoff = dropoff

def get_grav_accel(orbiter, body, dropoff=2):
    global grav_const
    
    grav_dir = (body.pos - orbiter.pos).normalized()
    grav_mag = grav_const * body.mass / ((body.pos - orbiter.pos).mag()**dropoff)

    return grav_dir * grav_mag

# this is an 8th order symplectic integrator
def propagateYoshida8(body, orbiters, end_time=120000, dt=1):
    # - - - CONSTANTS - - -
    w1 = 0.311790812418427e0
    w2 = -0.155946803821447e1
    w3 = -0.167896928259640e1
    w4 = 0.166335809963315e1
    w5 = -0.106458714789183e1
    w6 = 0.136934946416871e1
    w7 = 0.629030650210433e0
    w0 = 1.65899088454396 # (1 - 2 * (w1 + w2 + w3 + w4 + w5 + w6 + w7))

    ds = [w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7]

    # cs = [w7 / 2, (w7 + w6) / 2, (w6 + w5) / 2, (w5 + w4) / 2,
    #           (w4 + w3) / 2, (w3 + w2) / 2, (w2 + w1) / 2, (w1 + w0) / 2,
    #           (w1 + w0) / 2, (w2 + w1) / 2, (w3 + w2) / 2, (w4 + w3) / 2,
    #           (w5 + w4) / 2, (w6 + w5) / 2, (w7 + w6) / 2, w7 / 2]

    cs = [0.3145153251052165, 0.9991900571895715, 0.15238115813844, 0.29938547587066, -0.007805591481624963,
          -1.619218660405435, -0.6238386128980216, 0.9853908484811935, 0.9853908484811935, -0.6238386128980216,
          -1.619218660405435, -0.007805591481624963, 0.29938547587066, 0.15238115813844, 0.9991900571895715,
          0.3145153251052165]

    positions = []
    
    time = 0
    cycle = 0
    while time <= end_time:

        positions.append([])

        for orbiter in orbiters:
            positions[-1].append(orbiter.pos)

        for i in range(15):
            for orbiter in orbiters:
                orbiter.pos += orbiter.vel * cs[i] * dt
                grav_accel = get_grav_accel(orbiter, body, orbiter.dropoff)
                orbiter.vel += grav_accel * ds[i] * dt

        for orbiter in orbiters:
            orbiter.pos += orbiter.vel * cs[15] * dt
            
        cycle += 1
        time = cycle * dt

    return positions

# this is a 1st order symplectic integrator - see above for an 8th order one
def propagateSymplecticEuler(body, orbiters, end_time=12000, dt=1):

    positions = []
    
    time = 0
    cycle = 0
    while time <= end_time:

        positions.append([])

        for orbiter in orbiters:
            positions[-1].append(orbiter.pos)
            grav_accel = get_grav_accel(orbiter, body, orbiter.dropoff)
            orbiter.vel += grav_accel * dt
            orbiter.pos += orbiter.vel * dt

        cycle += 1
        time = cycle * dt

    return positions

def plot_positions(body, orbiters, positions):
    num_orbiters = len(positions[0])

    xs = []
    ys = []

    for timestamp in positions:
        xs.append([])
        ys.append([])
        for i in range(num_orbiters):
            xs[-1].append(timestamp[i].x)
            ys[-1].append(timestamp[i].y)

    xs_sanitized = []
    ys_sanitized = []
    for i in range(num_orbiters):
        xs_sanitized.append([])
        ys_sanitized.append([])

        for t in range(len(xs)):
            xs_sanitized[-1].append(xs[t][i])
            ys_sanitized[-1].append(ys[t][i])

    fig, ax = plt.subplots()
    
    for i in range(num_orbiters):
        ax.plot(xs_sanitized[i], ys_sanitized[i], label="Orbiter " + str(i) + ", dropoff " + str(orbiters[i].dropoff))

    body_circle = plt.Circle((0, 0), body.radius, color="r")
    ax.add_patch(body_circle)
    ax.axis("equal")
    plt.legend()
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.title("Orbit Plot")
    plt.grid()
    plt.show()

def main():
    orbiter1 = orbiter(vec2(2e11, 0), vec2(0, 120000))
    orbiter2 = orbiter(vec2(2e11, 0), vec2(0, 2500), 2.3)
    orbiter3 = orbiter(vec2(2e11, 0), vec2(0, 150e3), 1.98) 
    CygnusX1 = body(21.2 * 2e30, 21 * 7e8)

    orbiters = [orbiter1, orbiter2, orbiter3]

    positions = propagateYoshida8(CygnusX1, orbiters, 1e10, 2e5)

    plot_positions(CygnusX1, orbiters, positions)

main()
            
