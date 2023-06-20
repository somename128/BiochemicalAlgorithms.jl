using FFTW
using BiochemicalAlgorithms
using DSP
using FourierTools

include("extract_max.jl")

N = Int32(10)

A = zeros(Float64, N, N)
B = zeros(Float64, N, N)

for i in 3:10
    A[i] = i
end

# B[3,5,9] = 1

# C = ifft(fft(A).*fft(B))
#C_max = extract_max(C)
#D = irfft(rfft(A).*rfft(B), N)
#D_max = extract_max(D)
#E = DSP.conv(A,B)
#E_max = extract_max(E)
# max = extract_max(ccorr(A,B;centered=false))
fft([1,0,0])


