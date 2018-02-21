

class Polynomial:
    def __init__(self, n, k=None):
        self.degree = n
        if k is None:
            self.koef = [1 for i in range(n + 1)]
        else:
            if len(k) > n:
                k = k[:n + 1]
            elif len(k) < n:
                k.extend([0.0 for i in range(n - len(k))])
            self.koef = [x for x in k]
            if len(self.koef) < self.degree + 1:
                self.koef.extend([0 for i in range(self.degree - len(self.koef) + 1)])

    def __imul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            for i in range(len(self.koef)):
                self.koef[i] *= other
        elif isinstance(other, Polynomial):
            diff = other.degree + self.degree
            ext = [0.0 for i in range(diff + 1)]

            for i in range(other.degree + 1):
                for j in range(self.degree+1):
                    buff = self.koef[j] * other.koef[i]
                    if buff == 0:
                        continue
                    else:
                        ext[diff - (self.degree - j + other.degree - i)] += buff
            self.degree = diff
            self.koef = ext
        else:
            raise Exception
        return self

    def __mul__(self, other):
        ret = Polynomial(self.degree, self.koef)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(len(ret.koef)):
                ret.koef[i] *= other
        elif isinstance(other, Polynomial):
            diff = other.degree + ret.degree
            ext = [0.0 for i in range(diff + 1)]

            for i in range(other.degree + 1):
                for j in range(ret.degree + 1):
                    buff = ret.koef[j] * other.koef[i]
                    if buff == 0:
                        continue
                    else:
                        ext[diff - (ret.degree - j + other.degree - i)] += buff
            ret.degree = diff
            ret.koef = ext
        else:
            raise Exception
        return ret

    def __pow__(self, n):
        ret = Polynomial(self.degree, self.koef)
        for i in range(n - 1):
            ret *= self
        return ret

    def __iadd__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            self.koef[-1] += other
        elif isinstance(other, Polynomial):
            if other.degree > self.degree:
                diff = other.degree - self.degree
                self.degree += diff
                ext = [0.0 for i in range(diff)]
                ext.extend(self.koef)
                self.koef = ext

            i = self.degree
            while i != -1:
                self.koef[i] += other.koef[i]
                i -= 1
        else:
            raise Exception
        return self

    def __add__(self, other):
        ret = Polynomial(self.degree, self.koef)
        if isinstance(other, int) or isinstance(other, float):
            ret.koef[-1] += other
        elif isinstance(other, Polynomial):
            if other.degree > ret.degree:
                diff = other.degree - ret.degree
                ret.degree += diff
                ext = [0.0 for i in range(diff)]
                ext.extend(ret.koef)
                ret.koef = ext

                i = ret.degree
                while i != -1:
                    ret.koef[i] += other.koef[i]
                    i -= 1
            elif other.degree < ret.degree:
                buf = other.koef
                diff = ret.degree - other.degree
                ext = [0.0 for i in range(diff)]
                ext.extend(buf)
                i = ret.degree
                while i != -1:
                    ret.koef[i] += ext[i]
                    i -= 1
            else:
                i = ret.degree
                while i != -1:
                    ret.koef[i] += other.koef[i]
                    i -= 1
        else:
            raise Exception
        return ret

    def __isub__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            self.koef[-1] += other
        elif isinstance(other, Polynomial):
            if other.degree > self.degree:
                diff = other.degree - self.degree
                self.degree += diff
                ext = [0.0 for i in range(diff)]
                ext.extend(self.koef)
                self.koef = ext

            i = self.degree
            while i != -1:
                self.koef[i] -= other.koef[i]
                i -= 1
        else:
            raise Exception
        return self

    def __sub__(self, other):
        ret = Polynomial(self.degree, self.koef)
        if isinstance(other, int) or isinstance(other, float):
            ret.koef[-1] -= other
        elif isinstance(other, Polynomial):
            if other.degree > ret.degree:
                diff = other.degree - ret.degree
                ret.degree += diff
                ext = [0.0 for i in range(diff)]
                ext.extend(ret.koef)
                ret.koef = ext
            elif other.degree < ret.degree:
                buf = other.koef
                diff = ret.degree - other.degree
                ext = [0.0 for i in range(diff)]
                ext.extend(buf)
            else:
                ext = other.koef

            i = ret.degree
            while i != -1:
                ret.koef[i] -= ext[i]
                i -= 1
        else:
            raise Exception
        return ret

    def up_degree(self) :
        """...or multiplication polynomial by x"""
        self.degree += 1
        self.koef.append(0.0)

    def get(self, x):
        """Substitution x"""
        res = 0.0
        for i in range(self.degree + 1):
            res += (self.koef[i] * (x ** (self.degree - i)))
        return res

    def integral(self):
        ret = Polynomial(self.degree, self.koef)
        ret.koef = [k / (self.degree - ret.koef.index(k) + 1) if k else k
                    for k in ret.koef]
        ret.up_degree()
        return ret

    def definite_integral(self, down, up):
        i = self.integral()
        return i.get(up) - i.get(down)

    def integral_variable_up(self, down):
        i = self.integral()
        return i - i.get(down)

    def __str__(self):
        i = self.degree
        prnt_str = ""
        if not i:
            prnt_str += "{0}".format(self.koef[0])
            return prnt_str
        for k in self.koef:
            if k:
                if not i:
                        prnt_str += " + {0}".format(k)
                else:
                    if i != self.degree:
                        prnt_str += " + "
                    prnt_str += "{0}x^{1}".format(k, i)
            i -= 1
        return prnt_str


if __name__ == "__main__":
    """Tests"""
    # l = Polynomial(2, [3, 4, 9])
    # p = Polynomial(4, [2, 0, 0, 6, 1])
    # print(l, "\n", p)
    #
    # print(l * l)
    #
    # l = l * l
    #
    # print(l.integral_variable_up(1))

    y = Polynomial(3, [1/3])
    x = Polynomial(1, [1])

    print(y.koef, x.koef)

    x2 = x ** 2

    y2 = y ** 2

    print(x2 + y2)



