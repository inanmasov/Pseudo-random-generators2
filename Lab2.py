import sympy
from scipy.linalg import hankel

def gen_s(p, s):
    new_bit = 0
    for i in range(0, len(p)):
        new_bit += s[i] * p[i]
    new_bit %= 2
    out_bit = s[len(s) - 1]
    s = [new_bit] + s[:len(s) - 1]
    return out_bit, s


def write_to_file(result):
    count = 1
    with open('out.txt', 'w') as f:
        for item in result:
            f.write("%s " % item)
            if count == 60:
                f.write("\n")
                count = 1
            else:
                count += 1


def ReadSeq(Len):
    with open("out.txt", "r") as fin:
        text = fin.read().split()[0:Len]
    s = [int(x) for x in text]
    return s


def translate_to_list(value):
    return [int(i) for i in bin(value)[2:]]


def translate_from_list(data):
    val = 0
    for i in range(0, len(data) - 1):
        val += data[i]
        val *= 2
    val += data[-1]
    return val


def berlekamp(s, k):
    m = [1]
    l = [0]

    if len(s) % 2 == 0:
        N = 2 * k
    else:
        N = 2 * k + 1

    for t in range(0, N):
        m_t = translate_to_list(m[t])
        e_t = s[t]

        for i in range(1, l[t] + 1):
            if i >= len(m_t):
                break
            e_t ^= m_t[i] * s[t - i]

        if e_t == 1:
            lam = -1
            while lam <= t - 1:
                if l[lam + 1] >= l[t]:
                    break
                lam += 1
            m_lam = 1
            if lam >= 0:
                m_lam = m[lam]
            m_lam = translate_to_list(m_lam)
            deg = t - lam
            for j in range(0, deg):
                m_lam.insert(0, 0)
            for j in range(0, len(m_t)):
                if j >= len(m_lam):
                    m_lam.append(m_t[j])
                    continue
                m_lam[j] ^= m_t[j]
            m.append(translate_from_list(m_lam))
        elif e_t == 0:
            m.append(m[t])

        if e_t == 1 and l[t] <= t / 2:
            l.append(t + 1 - l[t])
        else:
            l.append(l[t])

    degree = l[-1]
    polynom = translate_to_list(m[-1])

    print("Минимальный многочлен последовательности:")
    for i in range(0, len(polynom) - 1):
        if polynom[i] != 0:
            print("x^{} + ".format(degree), end="")
        degree -= 1
    print("{}".format(polynom[-1]))


def hankel_det(s, k):
    firstcol = s[: k]
    lastrow = s[k - 1: 2 * k - 1]
    a = hankel(firstcol, lastrow)
    m = sympy.Matrix(a)
    return m.det() % 2, m


def eval_order(s, k):
    while True:
        for i in range(k, 0, -1):
            res, mat = hankel_det(s, i)
            if res != 0:
                print(f'k = {i}, определитель = {res}')
                sympy.pprint(mat)
                return i
            elif res == 0:
                print(f'k = {i}, определитель = {res}')
                if i == 1:
                    return i - 1


if __name__ == '__main__':
    p = 19 # 124657
    s = 7
    L = 64
    k = 32

    p = [int(a) for a in bin(p)[3:]]
    s = [int(bit) for bit in bin(s)[2:]]
    s = [0] * (len(p) - len(s)) + s

    L = int(L)
    output = []
    for i in range(L):
        new_s, s = gen_s(p, s)
        output.append(new_s)
    write_to_file(output)

    seq = ReadSeq(L)
    print("Длина последовательности = {}".format(len(seq)))
    order = eval_order(seq, k)
    print("Оценка порядока снизу = {}".format(order))

    subseq = seq[0: 2 * k]
    berlekamp(subseq, k)
