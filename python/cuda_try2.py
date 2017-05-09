import pycuda.autoinit
import math
import pycuda.driver as drv
import numpy
from pycuda.compiler import SourceModule

mod = SourceModule("""
__global__ void multiply_them(float *dest, float *a, float *b)
{
  const int i = threadIdx.x;
  const int bi = blockIdx.x;
  
  int destination_offset = blockDim.x * bi + i;
  
  // dest[i] = a[i] * b[i];
  dest[destination_offset] = (float)i;
}

__global__ void timesTwoKernel( float *d_Operand1, float *d_Result1, int arraySize ) {
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx > arraySize) return;
    d_Result1[idx] = d_Operand1[idx] * d_Operand1[idx];
} 

""")
array_size = 60000

multiply_them = mod.get_function("multiply_them")

a = numpy.random.randn(array_size).astype(numpy.float32)
b = numpy.random.randn(array_size).astype(numpy.float32)

dest = numpy.zeros_like(a)
block_size = 32
grid_size = int((array_size+block_size - 1)/block_size)
# (nx+threadsPerBlock.x-1)/threadsPerBlock.x
multiply_them(
        drv.Out(dest), drv.In(a), drv.In(b),
        block=(block_size, 1, 1), grid=(grid_size, 1))

print(dest[31])
print(dest[32])
print(dest[33])
s = numpy.int32(array_size)

timesTwoKrenel = mod.get_function('timesTwoKernel')
a[0]=1
a[1]=2
a[2]=3
s=numpy.int32(1)

timesTwoKrenel(
    drv.In(a),
    drv.Out(dest),
    s,
    block=(block_size, 1, 1),
    grid=(grid_size, 1)
    )

print(dest[0:3])
