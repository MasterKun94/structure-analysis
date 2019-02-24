# import itertools


class Factory:
    def __init__(self, type):
        self.type = type

    def __set__(self, instance, value):
        instance.__dict__[self.type] = value


class TryYYY:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        self.dd = self.plasss()

    def plasss(self):
        print('plas')
        return self.a + self.b + self.c

    # dd = plasss()

    def aa(self):
        print('aa')
        return self.dd + self.a

    def bb(self):
        print('bb')
        return self.dd + self.b

    def cc(self):
        print('cc')
        return self.dd + self.c


def fib(n):
    # if n < 2:
    #     return 1
    # else:
    #     return fib(n - 1) + fib(n - 2)
    f1 = f2 = 1
    for k in range(n):
        f1, f2 = f2, f2 + f1
    return f2


class Countable:
    counter = 0

    def __init__(self):
        Countable.counter += 1

    @classmethod
    def get_count(cls):
        return Countable.counter


def ttttt(a, b, c, **kwargs):
    print(a, b, c, kwargs)

    # yield from iter(args)


if __name__ == '__main__':

    a = ['1', '2', '3', '4']
    li = ", ".join([i for i in a])
    print(isinstance(li, str))
    for i in a:
        print(i)

    import itertools
    c = 0
    for i in itertools.count(0, 1):
        print(i)

