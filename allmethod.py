# euler 13.1
def euler(x0, y0, h, x):
    y = y0
    while x0 < x:
        y = y + h * (x0 + y + x0 * y)  # this is the function
        x0 += h
    return round(y, 4)


#  Range kutta 13.2
def f(x, y):
    fx = -2 * x * y ** 2
    return round(fx, 4)


x0 = int(input())
y0 = int(input())
x = int(input())
# second order
h = x - x0
k1 = h * f(x0, y0)
k2 = h * f(x0 + h, y0 + k1)
# k3 = h * f(x0 + 0.5 * (h), y0 + 0.5*(k2)) ...if in fourth order

y = y0 + 0.5 * (k1 + k2)
print(round(y, 4))


# Common for all of 12.*
def f(x):
    fx = 0
    order = len(coef)
    for i in range(order):
        fx += coef[i] * x ** (order - i - 1)
    return round(fx, 4)


coef = list(map(int, input().split()))
a, b = map(int, input().split())
n = int(input())
h = (b - a) / n
x = []
for i in range(n + 1):
    x.append(round(a + i * h, 4))  # Find all x terms
y = []
for i in range(n + 1):
    y.append(f(x[i]))  # Find all y terms
res = y[0] + y[-1]

# trapezoidal rule 12.1
# ...
for i in range(1, n):
    res += 2 * y[i]  # sum of first and last + twice the sum of the middle
res = res * h / 2
print(round(res, 4))

# simpson's 1/3 rule 12.2
# ...
for i in range(1, n):
    res += 2 * y[i] if i % 2 == 0 else 4 * y[i]
res = res * h / 3
print(round(res, 4))

# simpsons 3/8 rule 12.3
# ...
for i in range(1, n):
    res += 2 * y[i] if i % 3 == 0 else 3 * y[i]

res = res * 3 * h / 8
print(round(res, 4))


# Common for all 11.* except Newton-Raphson
def f(x):
    fx = 0
    order = len(coef)
    for i in range(order):
        fx += coef[i] * x ** (order - i - 1)
    return fx


# bisection method 11.1
# ...
def bisection(f, a, b, tol=1e-6):
    while (b - a) / 2.0 > tol:
        midpoint = (b - a) / 2
        if f(midpoint) == 0:
            return midpoint
        elif f(a) * f(midpoint) < 0:
            b = midpoint
        else:
            a = midpoint

    return (a + b) / 2.0


coef = list(map(int, input().split()))
a, b = map(float, input().split())
roundoff = int(input())
print(round(bisection(f, a, b), roundoff))


# regula falsi 11.2
# ...
def regulaFalsi(f, a, b, tol=1e-6):
    while abs(b - a) > tol:
        c = (a * f(b) - b * f(a)) / (f(b) - f(a))

        if abs(f(c)) < tol:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2.0


coef = list(map(int, input().split()))
a, b = map(int, input().split())
roundoff = int(input())
print(round(regulaFalsi(f, a, b), roundoff))


# newton-raphson 11.3
def f_x(x):
    fx = 0
    dfx = 0
    order = len(coef)
    for i in range(order):
        fx += coef[i] * x ** (order - i - 1)
    for i in range(order - 1):
        dfx += coef[i] * (order - i - 1) * x ** (order - i - 2)
    return fx, dfx


coef = list(map(int, input().split()))
a, b = map(float, input().split())
roundoff = int(input())
prev = (a + b) / 2
while True:
    f, df = f_x(prev)
    x = prev - f / df
    nxt = round(x, roundoff)
    if prev == nxt:
        break
    else:
        prev = nxt
print(nxt)


# secant method 11.4
# ...
def secant(a, b, tol=1e-6):
    while True:
        c = (a * f(b) - b * f(a)) / (f(b) - f(a))
        if abs(f(c)) < tol:
            return c
        a, b = b, c


coef = list(map(int, input().split()))
a, b = map(float, input().split())
roundoff = int(input())
print(round(secant(a, b), roundoff))

# mullers method 11.5
# ...
import math


def muller(f, x0, x1, x2, tol=1e-6, maxiter=100):
    for i in range(maxiter):
        h1, h2 = x1 - x0, x2 - x1
        i1, i2 = f(x1) - f(x0), f(x2) - f(x1)
        A = (i1 / h1 - i2 / h2) * 1 / (h1 + h2)
        B = A * i2 / h2
        rad = math.sqrt(B ** 2 - 4 * A * f(x2))

        if abs(B + rad) > abs(B - rad):
            den = B + rad
        else:
            den = B - rad
        dxr = 2 * f(x2) / den
        x3 = x2 - dxr
        if abs(dxr) < tol:
            return x3
        x0, x1, x2 = x1, x2, x3
    return None


a, b, c = map(int, input().split())
roundoff = int(input())
print(round(muller(f, a, b, c), roundoff))


# kruskal algorithm 10.1
class disjointSet:
    def __init__(self, vertices):
        self.parent = {v: v for v in range(vertices)}

    def find_root(self, v):
        if self.parent[v] != v:
            self.parent[v] = self.find_root(self.parent[v])
        return self.parent[v]

    def union(self, v1, v2):
        p1 = self.find_root(v1)
        p2 = self.find_root(v2)
        self.parent[p1] = p2


def kruskal(verts, edges):
    edges.sort()
    mst = []
    obj = disjointSet(verts)
    for edge in edges:
        w, s, d = edge
        p1 = obj.find_root(s)
        p2 = obj.find_root(d)
        if p1 != p2:
            mst.append(edge)
            obj.union(s, d)
    return mst


nv, ne = map(int, input().split())
verts = list(range(nv))
edges = [tuple(map(int, input().split())) for i in range(ne)]
mst = kruskal(verts, edges)
total = 0
for edge in mst:
    total += edge[0]
print(total)

# prims algorithm 10.2
import heapq


def prim(adjList, start):
    mst = []
    visited = set()
    heap = [(0, start)]

    while heap:
        w, d = heapq.heappop(heap)
        if d not in visited:
            visited.add(d)
            mst.append((w, d))
            for w, n in adjList[d]:
                if n not in visited:
                    heapq.heappush(heap, (w, n))
    return mst


nv, ne = map(int, input().split())
edges = list(map(int, input().split()))
adjlist = {v: [] for v in range(ne)}

for w, v1, v2 in edges:
    adjlist[v1].append((w, v2))
    adjlist[v2].append((w, v1))

start = 0
mst = prim(adjlist, start)
total = 0
for e in mst:
    total += e[0]
print(total)

# Common to all 9.*
nv, ne = map(int, input().split())
edges = [list(map(int, input().split())) for i in range(ne)]
adj_list = {i: [] for i in range(nv)}

for edge in edges:
    adj_list[edge[0]].append(edge[1])

start = int(input())


# BFS 9.1
# ...
def Bfs(graph, start):
    visited = set()
    queue = [start]
    visited.add(start)

    while queue:
        vertex = queue.pop(0)
        print(vertex, end=" ")

        for neighbor in graph[vertex]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)


Bfs(adj_list, start)


# DFS 9.2
def Dfs(graph, start, visited=None):
    if visited is None:
        visited = set()
    visited.add(start)
    print(start, end=" ")

    for neighbor in graph[start]:
        if neighbor not in visited:
            Dfs(graph, neighbor, visited)

    return visited


Dfs(adj_list, start)


# connected graphs 5.3
def dfs(vertex, visited):
    visited[vertex - 1] = True
    for neighbor in adjlist[vertex]:
        if neighbor in visited:
            dfs(neighbor, visited)


nv, ne = map(int, input().split())
visited = [False for i in range(nv)]
adjlist = {v + 1: [] for v in range(nv)}
edges = [tuple(map(int, input().split())) for i in range(ne)]

for v1, v2 in edges:
    adjlist[v1].append(v2)
    adjlist[v2].append(v1)

component = 0
for i in range(1, nv + 1):
    if visited[i - 1] == False:
        component += 1
        dfs(i, visited)
print(component)
