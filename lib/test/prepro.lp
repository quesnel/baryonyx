minimize

a + b + c + d + 1000.5 e + 1000 t + x + y + z + w + c1 + c2 + c3 + c4 + c5 + c6
 + d1 + d2 + d3 + d4 + d5 + d6 - 950 t - 75 t + 15 t + 25 t + 1000.45 - 998 e - 0.5 ZZ

subject to
  w = 0
  z + y + x = 3
- a - b - c <= -1
  a + e + d <= 2
  b + e + d <= 2
  a - c - b <= 0
  a + t = 0
  a + b + t = 1
  c1 + c2 + c3 >= 0
  c4 + c5 + c6 <= 3
  d1 + d2 + d3 >= 3
  d4 + d5 + d6 <= 0
  -ZZ <= 0

bound
t = 0

binaries
a b c d e w x y z c1 c2 c3 c4 c5 c6 d1 d2 d3 d4 d5 d6 ZZ

end
