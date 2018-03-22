srand(1234)

T = 10000000
a = 1+rand(T,1)
b = randn(T,1)
total = 0

println("Multiplication")
@time dot(a, b)
@time sum(a .* b)
@time sum([a[i]*b[i] for i in 1:T])
@time for i in 1:T; total += a[i]*b[i]; end

println("Summation")
egy = ones(T,1)
total = 0.0
@time sum(a)
@time dot(egy, a)
@time for i in 1:T; total += a[i]; end

println("Product")
@time prod(a)
@time exp(sum(log(a)))
