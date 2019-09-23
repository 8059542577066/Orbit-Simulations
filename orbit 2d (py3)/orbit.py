import math


class Orbit:

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
        self.alpha = self.tiltedAngle(r0, v0)
        self.ap = self.a * (1 + self.e) - self.R
        self.pe = self.a * (1 - self.e) - self.R

    def dotProduct(self, a, b):
        a1, a2 = a
        b1, b2 = b
        return a1 * b1 + a2 * b2

    def crossProduct(self, a, b):
        a1, a2 = a
        b1, b2 = b
        return a1 * b2 - a2 * b1

    def norm(self, v):
        x, y = v
        return math.sqrt(x ** 2 + y ** 2)

    def trueAnomaly(self):
        return math.acos((self.H ** 2 / (self.mu * self.r0) - 1) / self.e)

    def tiltedAngle(self, r0, v0):
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
        print("r = " + str(self.H ** 2) + " / (1000 * " + str(self.mu) +\
              " * (1 + " + str(self.e) + " * cos(\\theta - " +\
              str(self.alpha) + ")))")


orbit = Orbit(3.986004418E+14, 6371000, (6000000, -8000000), (4000, 5000))
