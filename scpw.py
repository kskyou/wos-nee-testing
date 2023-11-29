# Slowest Closest Points in the West
import math
import numpy as np

PI = 3.1415926

def dist_segment(point, segment): 
    # https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    x1,y1 = segment[0]
    x2,y2 = segment[1]
    x3,y3 = point
    px = x2-x1
    py = y2-y1
    norm = px*px + py*py
    u =  ((x3 - x1) * px + (y3 - y1) * py) / norm
    if u > 1:
        u = 1
    elif u < 0:
        u = 0
    x = x1 + u * px
    y = y1 + u * py
    dx = x - x3
    dy = y - y3
    dist = (dx*dx + dy*dy)**(1/2)
    return ((x,y), dist)


def closest_point(point, boundary):
    bestpoint = (0, 0)
    bestdist = 1000000
    for segment in boundary:
        newpoint, newdist = dist_segment(point, segment)
        if (newdist < bestdist):
            bestdist = newdist
            bestpoint = newpoint
    return (bestpoint, bestdist)


def create_boundary(L):
    assert(len(L) > 2)
    boundary = [(L[-1], L[0])]
    for i in range(len(L)-1):
        boundary.append((L[i],L[i+1]))
    return boundary


# Ray boundary intersection. Currently not needed. 
'''

def magnitude(vector):
   return np.sqrt(np.dot(np.array(vector),np.array(vector)))

def norm(vector):
   return np.array(vector)/magnitude(np.array(vector))

def ray_segment(rayOrigin, rayDirection, segment):
    # https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin
    point1, point2 = segment
    point1 = list(point1)
    point2 = list(point2)
    rayOrigin = list(rayOrigin)
    rayDirection = list(rayDirection)

    rayOrigin = np.array(rayOrigin, dtype=float)
    rayDirection = np.array(norm(rayDirection), dtype=float)
    point1 = np.array(point1, dtype=float)
    point2 = np.array(point2, dtype=float)
    v1 = rayOrigin - point1
    v2 = point2 - point1
    v3 = np.array([-rayDirection[1], rayDirection[0]])
    t1 = np.cross(v2, v1) / np.dot(v2, v3)
    t2 = np.dot(v1, v3) / np.dot(v2, v3)
    if t1 >= 0.0 and t2 >= 0.0 and t2 <= 1.0:
        intersection = (rayOrigin[0] + t1 * rayDirection[0], rayOrigin[1] + t1 * rayDirection[1])
        return (intersection, t1)
    return None

def ray_boundary(point, ray, boundary):
    best = None
    for segment in boundary:
        new = ray_segment(point, ray, segment)
        if (new != None):
            if (best == None or new[1] < best[1]):
                best = new
    return best

def near_boundary(point, distance, boundary):
    L = []
    for i in range(len(boundary)):
        segment = boundary[i]
        _, newdist = dist_segment(point, segment)
        if (newdist < distance):
            L.append(i)

    return L
'''
