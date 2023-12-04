import math
import numpy as np
import scpw
import sampling

EPSILON = 0.005
MAX_LENGTH = 1000
nWalks = 1000
RUNS = 200

boundary_square = scpw.create_boundary([(0,0),(1,0),(1,1),(0,1)]) # square
boundary_keyhole = scpw.create_boundary([(0,-0.5),(1,-0.5),(1,-1),(2,-1),(2,1),(1,1),(1,0.5),(0,0.5)]) #keyhole

BIG = 10000000
boundary_half = scpw.create_boundary([(-BIG,0),(BIG,0),(BIG,BIG),(-BIG,BIG)]) # half plane

f1 = lambda p : 1 # uniform
f2 = lambda p : -math.log((p[0]+0.1) ** 2 + (p[1]) ** 2) # Point source at (-0.1, 0)
f3 = lambda p : math.exp(3 * p[0]) * math.sin(3 * p[1]) # holomorphic
f4 = lambda p : 0 if p[0] < 1 else math.sin(p[0] * 20) ** 2 # high frequency

point_middle = (0.5, 0.5)
point_edge = (1.1 * EPSILON, 0.2)


def walk_tiny(problem):
    boundary = problem[0]
    curr = problem[1]
    f = problem[2]

    curr_length = 0
    p, r = scpw.closest_point(curr, boundary)

    total = 0
    multiplier = 1
    while (curr_length < MAX_LENGTH):
        curr_length += 1

        pdirection = sampling.normalize((p[0] - curr[0], p[1] - curr[1]))  # Direction to closest point
        angle = 2 * math.sin((EPSILON/2) / r) # Ball-ball intersection
        angle = min(angle, 0.1)

        kernel = 1 / (2 * math.pi * r)
        pc = 1 / (2 * angle * r) # Probability of cap
        pb = 1 / (2 * (math.pi - angle) * r) # Probability of ball minus cap

        qdirection = sampling.sample_cap(pdirection, angle) # Sampled direction
        q = (curr[0] + qdirection[0] * r, curr[1] + qdirection[1] * r) # Sampled point
        qp, qr = scpw.closest_point(q, boundary)

        #assert(qr < EPSILON)
        total += f(qp) / pc * kernel * multiplier

        direction = sampling.sample_cap((-pdirection[0],-pdirection[1]), math.pi - angle ) #Sampled direction not in cap
        curr = (curr[0] + direction[0] * r, curr[1] + direction[1] * r)
        p, r = scpw.closest_point(curr, boundary)

        #assert (sampling.dist(pdirection, direction) >= 2 * math.sin(angle/2))
        multiplier *= kernel / (pb)

        if (r < EPSILON): # If in epsilon shell
            return (total + f(p) * multiplier, curr_length)
    
    print("Exceeded max length")
    return 0


'''
def walk_mises_mis(problem):
    boundary = problem[0]
    curr = problem[1]
    f = problem[2]

    curr_length = 0
    p, r = scpw.closest_point(curr, boundary)

    total = 0
    multiplier = 1
    while (curr_length < MAX_LENGTH):

        pdirection = sampling.normalize((p[0] - curr[0], p[1] - curr[1]))  # Direction to closest point
        #angle = math.cos((r - EPSILON) / r) # Half width of cap
        angle = 1

        kernel = 1 / (2 * math.pi * r)
        circlelength = 2 * math.pi * r # Circumference of circle

        qdirection = sampling.sample_mises(pdirection, angle) # Sampled direction
        q = (curr[0] + qdirection[0] * r, curr[1] + qdirection[1] * r) # Sampled point

        qp, qr = scpw.closest_point(q, boundary)
        if (qr < EPSILON): # If in epsilon shell
            total += f(qp) / (sampling.pdf_mises(pdirection, angle, r, qdirection) + 1 / circlelength) * kernel * multiplier

        oldr = r
        direction = sampling.sample_circle()
        curr = (curr[0] + direction[0] * r, curr[1] + direction[1] * r)
        p, r = scpw.closest_point(curr, boundary)

        if (r < EPSILON): # If in epsilon shell
            multiplier *= kernel / (sampling.pdf_mises(pdirection, angle, oldr, direction) + 1 / circlelength)

        if (r < EPSILON): # If in epsilon shell
            return total + f(p) * multiplier
'''

    
def walk_cap_mis(problem):
    boundary = problem[0]
    curr = problem[1]
    f = problem[2]

    curr_length = 0
    p, r = scpw.closest_point(curr, boundary)

    total = 0
    multiplier = 1
    while (curr_length < MAX_LENGTH):
        curr_length += 1

        pdirection = sampling.normalize((p[0] - curr[0], p[1] - curr[1]))  # Direction to closest point
        angle = math.cos((r - EPSILON) / r) # Half width of cap

        kernel = 1 / (2 * math.pi * r)
        caplength = 2 * angle * r # Length of cap
        circlelength = 2 * math.pi * r # Circumference of circle

        qdirection = sampling.sample_cap(pdirection, angle) # Sampled direction
        q = (curr[0] + qdirection[0] * r, curr[1] + qdirection[1] * r) # Sampled point
        qp, qr = scpw.closest_point(q, boundary)
        if (qr < EPSILON): # If in epsilon shell
            total += f(qp) / (1 / caplength + 1 / circlelength) * kernel * multiplier

        oldr = r
        direction = sampling.sample_circle()
        curr = (curr[0] + direction[0] * r, curr[1] + direction[1] * r)
        p, r = scpw.closest_point(curr, boundary)

        if (sampling.dist(pdirection, direction) < 2 * math.sin(angle/2) and r < EPSILON): # If in cap and epsilon shell
            multiplier *= kernel / (1 / caplength + 1 / circlelength)

        if (r < EPSILON): # If in epsilon shell
            return (total + f(p) * multiplier, curr_length)
    
    print("Exceeded max length")
    return 0


def walk_cap(problem):
    boundary = problem[0]
    curr = problem[1]
    f = problem[2]

    curr_length = 0
    p, r = scpw.closest_point(curr, boundary)

    total = 0
    while (curr_length < MAX_LENGTH):
        curr_length += 1

        pdirection = sampling.normalize((p[0] - curr[0], p[1] - curr[1]))  # Direction to closest point
        angle = math.cos((r - EPSILON) / r) # Half width of cap

        kernel = 1 / (2 * math.pi * r)
        caplength = 2 * angle * r

        qdirection = sampling.sample_cap(pdirection, angle) # Sampled direction
        q = (curr[0] + qdirection[0] * r, curr[1] + qdirection[1] * r) # Sampled point
        qp, qr = scpw.closest_point(q, boundary)
        if (qr < EPSILON): # If in epsilon shell
            total += f(qp) * caplength * kernel

        oldr = r
        direction = sampling.sample_circle()
        curr = (curr[0] + direction[0] * r, curr[1] + direction[1] * r)
        p, r = scpw.closest_point(curr, boundary)
        if (r < EPSILON): # If in epsilon shell
            if (sampling.dist(pdirection, direction) < 2 * math.sin(angle/2)): # If in cap
                return (total, curr_length)
            else:
                return (total + f(p), curr_length)
    
    print("Exceeded max length")
    return 0


def walk_naive(problem):
    boundary = problem[0]
    curr = problem[1]
    f = problem[2]

    curr_length = 0
    p, r = scpw.closest_point(curr, boundary)

    while (curr_length < MAX_LENGTH):
        curr_length += 1
        
        direction = sampling.sample_circle()
        curr = (curr[0] + direction[0] * r, curr[1] + direction[1] * r)
        p, r = scpw.closest_point(curr, boundary)

        if (r < EPSILON):
            return (f(p), curr_length)

    print("Exceeded max length")
    return 0


def estimate_solution(walk, nWalks, problem):
    total = 0
    total_l = 0
    for w in range(nWalks):
        res, length = walk(problem)
        total += res
        total_l += length
    return (total / nWalks, total_l / nWalks)


if __name__ == "__main__":

    ss = [1, 2, 4, 8, 16, 32, 64]

    for s in ss:
        print("---------- s = " + str(s) + "----------")
        f5 = lambda p : 10 * s if (p[0] > 0.5 - 0.5 / s  and p[0] < 0.5 + 0.5 / s) else 0 # point 
        problem = (boundary_square, point_middle, f5)

        #print("true (if f harmonic)", problem[2](problem[1]))
        walks = [walk_naive, walk_tiny, walk_cap, walk_cap_mis]
        for walk in walks:
            print(walk.__name__)

            lengths = []
            for i in range(10):
                res, length = estimate_solution(walk, nWalks, problem)
                lengths.append(length)
            lengths = np.array(lengths)
            #print("avg len", np.average(lengths))
            nWalks_adj = nWalks * 10 / np.average(lengths)

            results = []
            for i in range(RUNS):
                res, length = estimate_solution(walk, int(nWalks_adj), problem)
                results.append(res)
            results = np.array(results)
            print("mean", np.average(results))
            print("std", np.std(results))

