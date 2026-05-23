#include <iostream>
#include <TApplication.h>
#include <TH1F.h>
#include <TCanvas.h>

// CUDA kernel for adding arrays
__global__ void addArrays(const float* a, const float* b, float* result, int size) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < size) {
        result[idx] = a[idx] + b[idx];
    }
}

int main(int argc, char** argv) {
    const int size = 100;
    float h_a[size], h_b[size], h_result[size];
    float *d_a, *d_b, *d_result;

    // Initialize host arrays
    for (int i = 0; i < size; i++) {
        h_a[i] = i * 0.5;
        h_b[i] = i * 1.0;
    }

    // Allocate device memory and copy data
    cudaMalloc((void**)&d_a, size * sizeof(float));
    cudaMalloc((void**)&d_b, size * sizeof(float));
    cudaMalloc((void**)&d_result, size * sizeof(float));
    cudaMemcpy(d_a, h_a, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, size * sizeof(float), cudaMemcpyHostToDevice);

    // Launch CUDA kernel
    addArrays<<<(size + 255) / 256, 256>>>(d_a, d_b, d_result, size);

    // Copy result back to host
    cudaMemcpy(h_result, d_result, size * sizeof(float), cudaMemcpyDeviceToHost);

    // Initialize ROOT application and plot result
    TApplication app("app", &argc, argv);
    TCanvas canvas("canvas", "CUDA and ROOT Example", 800, 600);
    TH1F hist("hist", "Array Addition Result", size, 0, size);
    for (int i = 0; i < size; i++) {
        hist.SetBinContent(i+1, h_result[i]);
    }
    hist.Draw();
    canvas.Update();
    app.Run();

    // Free device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_result);

    return 0;
}
