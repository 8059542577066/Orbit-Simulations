import math


class Orbit2D:

    def __init__(self, mu, R, r0, v0):
        self.mu = mu
        self.R = R
        # these are saved to draw tangent testing lines.
        self.x0, self.y0 = r0
        self.x0, self.y0 = self.x0 / 1000, self.y0 / 1000
        self.vx, self.vy = v0
        # end of tangent line arguments
        self.r0 = self.norm(r0)
        self.v0 = self.norm(v0)
        self.H = abs(self.crossProduct(r0, v0))
        self.a = mu * self.r0 / (2 * mu - self.r0 * self.v0 ** 2)
        self.e = math.sqrt(1 - self.H ** 2 / (mu * self.a))
        self.peTheta = self.periapsisAngle(r0, v0)
        self.ap = self.a * (1 + self.e) - self.R
        self.pe = self.a * (1 - self.e) - self.R

    def norm(self, v):
        x, y = v
        return math.sqrt(x ** 2 + y ** 2)

    def dotProduct(self, a, b):
        a1, a2 = a
        b1, b2 = b
        return a1 * b1 + a2 * b2

    def crossProduct(self, a, b):
        a1, a2 = a
        b1, b2 = b
        return a1 * b2 - a2 * b1

    def trueAnomaly(self):
        return math.acos((self.H ** 2 / (self.mu * self.r0) - 1) / self.e)

    def periapsisAngle(self, r0, v0):
        x0, y0 = v0
        vx, vy = v0
        rad = math.atan2(self.y0, self.x0)
        if self.dotProduct(r0, v0) > 0:
            return rad - self.trueAnomaly()
        else:
            return rad + self.trueAnomaly()

    def printTangent(self):
        print("y = " + str(self.y0) + " + (" + str(self.vy) + " / " +\
              str(self.vx) + ") * (x - " + str(self.x0) + ")")

    def printDesmos(self):
        print("r = " + "%f" % self.H ** 2 + " / (1000 * " + str(self.mu) +\
              " * (1 + " + str(self.e) + " * cos(\\theta - " +\
              str(self.peTheta) + ")))")


class Orbit3D:

    def __init__(self, mu, R, r0, v0):
        self.phi, self.theta = self.getAngles(self.crossProduct(r0, v0))
        r02d, v02d = self.antiRotate([r0, v0], self.phi, self.theta)
        r02d, v02d = r02d[:-1], v02d[:-1]
        orbit = Orbit2D(mu, R, r02d, v02d)
        # testing if 2D orbit works
        self.H, self.mu, self.e = orbit.H, orbit.mu, orbit.e
        # I changed below function to use differential integration
        #points = self.getTrajectory(r0, v0)
        points = self.getOrbit(orbit.peTheta)
        points += self.getTangent(r0, v0)
        self.save(points)
        # end of testing

    def crossProduct(self, a, b):
        a1, a2, a3 = a
        b1, b2, b3 = b
        return a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1

    def getAngles(self, vector):
        x, y, z = vector
        if (x, y) == (0, 0):
            return 0, 0
        else:
            rho = math.sqrt(x ** 2 + y ** 2 + z ** 2)
            phi = math.acos(z / rho)
            theta = math.atan2(y, x)
        return phi, theta

    def rotateY(self, points, rad):
        rotated = []
        for x, y, z in points:
            rotated.append((x * math.cos(rad) + z * math.sin(rad),\
                            y,\
                            z * math.cos(rad) - x * math.sin(rad)))
        return rotated

    def rotateZ(self, points, rad):
        rotated = []
        for x, y, z in points:
            rotated.append((x * math.cos(rad) - y * math.sin(rad),\
                            x * math.sin(rad) + y * math.cos(rad),\
                            z))
        return rotated

    def rotate(self, points, phi, theta):
        points = self.rotateY(points, phi)
        points = self.rotateZ(points, theta)
        return points

    def antiRotate(self, points, phi, theta):
        points = self.rotateZ(points, -theta)
        points = self.rotateY(points, -phi)
        return points

    def getPoint(self, rad):
        r = self.H ** 2 / (self.mu * (1 + self.e * math.cos(rad)))
        return r * math.cos(rad), r * math.sin(rad), 0

    def getOrbit(self, peTheta):
        points = []
        for i in range(2 ** 16):
            rad = 2 * math.pi * i / 2 ** 16
            points.append(self.getPoint(rad))
        # rotation testing
        points = self.rotateZ(points, peTheta)
        points = self.rotate(points, self.phi, self.theta)
        # printed in mega meter unit
        points = [(x / 10 ** 6, y / 10 ** 6, z / 10 ** 6) for\
                  (x, y, z) in points]
        return points

    # Testing if velocity vector is aligned
    def getTangent(self, r0, v0):
        x0, y0, z0 = r0
        x0, y0, z0 = x0 / 10 ** 6, y0 / 10 ** 6, z0 / 10 ** 6
        vx, vy, vz = v0
        vx, vy, vz = vx / 10 ** 3, vy / 10 ** 3, vz / 10 ** 3
        points = []
        for i in range(-5000, 15000):
            points.append((i / 20000 * vx + x0,\
                           i / 20000 * vy + y0,\
                           i / 20000 * vz + z0))
        return points

    # Integral numerical analysis of orbit trajectory
    def getForceField(self, r):
        x, y, z = r
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        return -self.mu / r ** 3 * x,\
               -self.mu / r ** 3 * y,\
               -self.mu / r ** 3 * z

    def getNextVelocity(self, v, f, dt):
        vx, vy, vz = v
        fx, fy, fz = f
        return vx + fx * dt, vy + fy * dt, vz + fz * dt

    def getNextPos(self, r, v, dt):
        x, y, z = r
        vx, vy, vz = v
        return x + vx * dt, y + vy * dt, z + vz * dt

    def getTrajectory(self, r0, v0):
        points = []
        points.append(r0)
        fx, fy, fz = self.getForceField(r0)
        dt = 1 / 2 ** 4  # differential to integrate
        v = v0
        r = r0
        for i in range(2 * 10 ** 5):
            v = self.getNextVelocity(v, self.getForceField(r), dt)
            r = self.getNextPos(r, v, dt)
            points.append(r)
        points = [(x / 10 ** 6, y / 10 ** 6, z / 10 ** 6) for\
                  (x, y, z) in points]
        return points
    # End of differential integration method

    def save(self, points):
        with open("output.txt", "w") as f:
            f.write("import bpy\n\n\n")
            f.write("vertices = []\n\n")
            for point in points:
                f.write("vertices.append(" + str(point) + ")\n")
            f.write("\n\nmesh = bpy.data.meshes.new(\"orbit\")\n")
            f.write("object = bpy.data.objects.new(\"orbit\", mesh)\n\n")
            f.write("object.location = bpy.context.scene.cursor_location\n")
            f.write("bpy.context.scene.objects.link(object)\n\n")
            f.write("mesh.from_pydata(vertices, [], [])\n")
            f.write("mesh.update()\n")


# orbit = Orbit2D(3.986004418E+14, 6371000, (6000000, -8000000), (4000, 5000))
orbit = Orbit3D(3.986004418E+14, 6371000,\
                (-4471000, -3000000, 2000000), (-2000, 7800, 3000))
