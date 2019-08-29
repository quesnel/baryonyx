minimize

a + b + c + d + 1000.5 ve + 1000 t + x + y + z + w + c1 + c2 + c3 + c4 + c5 + c6
 + d1 + d2 + d3 + d4 + d5 + d6 - 950 t - 75 t + 15 t + 25 t - 998 ve - 0.5 ZZ

subject to
ct00:  w = 0
ct01:  z + y + x = 3
ct02:- a - b - c <= -1
ct03:  a + ve + d <= 2
ct04:  b + ve + d <= 2
ct05:  a - c - b <= 0
ct06:  a + t = 0
ct07:  a + b + t = 1
ct08:  c1 + c2 + c3 >= 0
ct09:  c4 + c5 + c6 <= 3
ct10:  d1 + d2 + d3 >= 3
ct12:  d4 + d5 + d6 <= 0
ct13:  -ZZ <= 0

bound
t = 0

binaries
a b c d ve w x y z c1 c2 c3 c4 c5 c6 d1 d2 d3 d4 d5 d6 ZZ t

end

