// This is an implementation of the radix-2-DIT FFT algorithm.
// As of now the function is very much unsafe as the algorithm only is valid for arrays of size 2^N.
// The function has been tested and works with no std.
// X = fft(x)
// Having size 1 as bottom layer might not be optimal.

// x is an array of size N. x is the signal on which the FFT is performed.
// X is an array of size 2*N, because of complex numbers. X = [r0,c0,r1,c1...rN-1,cN-1]. X is the result of the fourier transform.
// twiddle_factor is an array of size N containing complex exponentials twiddle_factor = [cos(-2*pi*0/N), sin(-2*pi*0/N), cos(-2*pi*1/N), sin(-2*pi*1/N)...cos(-2*pi*(N/2-1)/N), sin(-2*pi*(N/2-1)/N)]
// In general the array sizes should be seen as minimum array sizes. If N = 64, an x with 0 elements will still work, an x with 55 elements however, will not work.
// N is the length of the array, has to be power of 2.
// s should be 1 when called by the user. It lets the same twiddle factor array to be used all the way down, taking every s:th element. It is also used when doing the odd side recursion.
unsafe fn ditfft2( x: &[f32], X: &mut [f32], twiddle_factor: &[f32], N: u32, s: u32 ) {
    let PI = 3.14159265; // There is a way to do this better. Constants are a thing.
    if N == 1 {
        X[0]=x[0];
    }
    else {
        ditfft2( &x, &mut X[..N as usize], &twiddle_factor, N/2, 2*s); //should probably bitshift N and s. even elements.
        ditfft2( &x[s as usize..], &mut X[N as usize..], &twiddle_factor, N/2, 2*s); // slice. nice. Bitshifting n and s still applies. odd elements.

        for k in (0..N).filter(|&k| k % 2 == 0) {
            let even_real = X[k as usize];
            let even_complex = X[(k+1) as usize];
            let odd_real = X[(N+k) as usize];
            let odd_complex = X[(N+k+1) as usize];
            X[k as usize] = even_real+twiddle_factor[(s*k) as usize]*odd_real-twiddle_factor[(s*k+1) as usize]*odd_complex; // By doin some clever premultiplication 4 multiply/iteration could be saved.
            X[(k+1) as usize] = even_complex +twiddle_factor[(s*k+1) as usize]*odd_real+twiddle_factor[(s*k) as usize]*odd_complex;
            X[(k+N) as usize] = even_real - twiddle_factor[(s*k) as usize]*odd_real + twiddle_factor[(s*k+1) as usize]*odd_complex;
            X[(k+N+1) as usize] = even_complex - twiddle_factor[(s*k) as usize]*odd_complex - twiddle_factor[(s*k+1) as usize]*odd_real;
        }
    }

}
