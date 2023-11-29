import math
import random
from scipy.stats import vonmises

K = 2

def randf():
    return random.uniform(0, 1)

def sample_circle():
    theta = randf() * 2 * math.pi
    return (math.cos(theta), math.sin(theta))

def sample_cap(direction, angle):
    phi = math.atan2(direction[1], direction[0])
    theta = (2 * randf() - 1) * angle
    return (math.cos(theta+phi), math.sin(theta+phi))

def sample_mises(direction, angle):
    angle = max(angle, 0.01)
    loc = math.atan2(direction[1], direction[0])
    kappa = K / (angle * angle)
    phi = (vonmises(loc=loc, kappa=kappa).rvs(1))[0]
    return (math.cos(phi), math.sin(phi))

def pdf_mises(direction, angle, r, pdirection):
    angle = max(angle, 0.01)
    loc = math.atan2(direction[1], direction[0])
    kappa = K / (angle * angle)
    ploc = math.atan2(pdirection[1], pdirection[0])
    return vonmises.pdf(loc=loc, kappa=kappa, x=ploc) / r 

def normalize(p):
    return (p[0] / math.sqrt(p[0] ** 2 + p[1] ** 2), p[1] / math.sqrt(p[0] ** 2 + p[1] ** 2))

def dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

'''
def sample_segment(segment):
    p1, p2 = segment
    length = dist(p1, p2) * randf()
    return (p1[0] + length * (p2[0] - p1[0]), p1[1] + length * (p2[1] - p1[1]))
'''

'''
total = 0
angle = 0.1
for i in range(10000):
    direction = sample_mises((1, 0), angle)
    total += 1 / pdf_mises((1, 0), angle, 1, direction)
print(total / 10000 )
'''

