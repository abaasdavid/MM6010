from sympy import Rational, factorint

class Qp:
    
    def __init__(self, m, p, prec=10):
        self.p = p
        self.prec = prec
        m = Rational(m, 1)
        self.v = Qp.valuation(m, p)  
        self.digits = self._digits(m)
        self.m = self._new_m()
    
    @staticmethod
    def valuation(m, p):  
        """Beräknar valueringen, returnerar None om talet 
           är 0 modulo precisionen"""
        num, den = m.p, m.q
        if num == 0:
            return None
        a = factorint(num, limit=p).get(p, 0)
        b = factorint(den, limit=p).get(p, 0)
        return a - b

    def _digits(self, m):
        digits = {}
        if not self:
            return {self.prec-1: 0}
        m0 = (Rational(self.p, 1)**(-self.v))*m
        k_n = m0.p
        b_prime = m0.q
        for n in range(self.v, self.prec):
            a_n = (k_n * pow(b_prime, -1, self.p)) % self.p
            k_n = (k_n - b_prime * a_n) // self.p
            digits[n] = a_n
        return digits
    
    def _new_m(self): 
        """Skapar en ny representation av det p-adiska talet
           med minsta möjliga potenser av p"""
        new_m = 0
        for power, digit in self.digits.items():
            if power < self.prec: 
                new_m += digit * (Rational(self.p, 1) ** power)
        return Rational(new_m, 1)
    
    def __repr__(self):
        return (f"Qp(Rational({self.m.p}, {self.m.q}), "
                f"p = {self.p}, prec={self.prec})")
    
    def _latex(self):
        if self.m == 0:
            return f"O\\left({self.p}^{{ {self.prec} }}\\right)"
        p_adic_terms = [ f"{digit} \\cdot {self.p}^{{ {power} }}" 
                        if power != 0 else str(digit)
                        for power, digit in sorted(self.digits.items())
                        if digit != 0 ]
        highest_power = max(self.digits.keys(), default=-1)
        return (f"{' + '.join(p_adic_terms)} "
                f"+ O\\left({self.p}^{{ {highest_power + 1} }}\\right)")

    def _repr_latex_(self):
        return f"$ {self._latex()} $" 

    def __str__(self):
        compact_str = ['...']
        separator = '_' if self.p >= 11 else ''
        highest_power = max(self.digits.keys(), default=0)
        for power in range(highest_power, -1, -1):
            compact_str.append(str(self.digits.get(power, 0)))
            if power != 0:
                compact_str.append(separator)
        lowest_power = min(self.digits.keys(), default=-1)
        if lowest_power < 0:
            compact_str.append(',')
            for power in range(-1, lowest_power - 1, -1):
                digit = self.digits.get(power)
                if digit is None:
                    compact_str.append('.')
                else:
                    compact_str.append(str(digit))
                compact_str.append(separator)
        return ''.join(compact_str).rstrip(separator)

    def __add__(self, other):
        if not isinstance(other, Qp):
            return Qp(self.m + Rational(other,1), self.p, prec=self.prec)
        else:
            if not self.p == other.p:
                raise ValueError("Addition av p-adiska tal med olika baser.")
            min_prec = min(self.prec, other.prec)
            new_m = self.m + other.m
            return Qp(new_m, self.p, prec=min_prec)
        
    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Qp(-self.m, self.p, prec=self.prec)
    
    def __pos__(self):
        return self
    
    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)
    
    def __mul__(self, other):
        if not isinstance(other, Qp):
            other = Rational(other, 1)
            if other == 0: 0
            v = Qp.valuation(other, self.p)
            return Qp(self.m * other, self.p, prec=self.prec + v)
        if not self.p == other.p:
            raise ValueError("Multiplikation av p-adiska tal med olika baser.")
        if not self or not other:                                                   
            max_prec = max(self.prec, other.prec)
            return Qp(0, self.p, max_prec)
        min_prec = min(self.v + other.prec, other.v + self.prec)
        return Qp(self.m * other.m, self.p, prec=min_prec) 

    def __rmul__(self, other):
        return self * other
            
    def __invert__(self):
        if not self:
            raise ZeroDivisionError 
        min_prec = -self.v + self.prec
        return Qp(1 / self.m, self.p, prec=min_prec)
    
    def __pow__(self, n):
        if n < 0 and not self:
            raise ZeroDivisionError("Kan inte upphöja 0 till ett negativt tal.")
        if n == 0:
            return Qp(1, self.p, prec=self.prec)
        if not self:
            return Qp(0, self.p, prec=self.prec ** n)
        new_m = self.m ** n
        new_v = self.v * n
        new_prec = self.prec + (new_v - self.v)
        return Qp(new_m, self.p, prec=new_prec)
    
    def __truediv__(self, other):
        if not isinstance(other, Qp):
            if other == 0: raise ZeroDivisionError
            return self * Rational(1,other) # sympy tycker ej om Rational(1,0)..
        return self * ~other

    def __rtruediv__(self, other):
        return ~self * other

    def __eq__(self, other):
        if not isinstance(other, Qp):
            other = Qp(other, self.p, prec=self.prec)
        if self.p != other.p:
            return False
        min_prec = min(self.prec, other.prec)
        for i in range(min_prec):
            if self.digits.get(i, 0) != other.digits.get(i, 0):
                return False
        return True
    
    def __bool__(self):
        return self.v is not None
    
    def __abs__(self):
        return 0 if not self else Rational(self.p, 1)**(-self.v)