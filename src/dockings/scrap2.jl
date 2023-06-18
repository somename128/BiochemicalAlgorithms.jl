using FFTW
using BiochemicalAlgorithms
using DSP
using FourierTools

include("extract_max.jl")

N = 10

A = zeros(Float32, N, N, N)
B = zeros(Float32, N, N, N)

for i in 1:N*N*N
    A[i] = i
end

B[3,5,9] = 1

C = ifft(fft(A).*fft(B))
#C_max = extract_max(C)
#D = irfft(rfft(A).*rfft(B), N)
#D_max = extract_max(D)
#E = DSP.conv(A,B)
#E_max = extract_max(E)
ccorr(A,B;centered=false)

