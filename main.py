v, e = map(int, input().split(', '))
li = [list(map(int, input().split(', '))) for i in range(e)]


def adjmat(v, e):
    mat = [[0] * v for i in range(v)]
    for edge in li:
        x, y = edge
        mat[x][y] = 1
        mat[y][x] = 1
    return mat


def adjlist(v, e):
    matlis = {}
    for i in range(v):
        matlis[i] = []

    for edge in li:
        x, y = edge
        matlis[x].append(y)
    for ver, neigh in matlis.items():
        print(f"{ver} : {neigh}")


print("Adjacency Matrix")
for i in (adjmat(v, e)):
    print(*i, sep=", ")
print("Adjacency List")
(adjlist(v, e))
