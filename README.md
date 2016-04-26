# NBody-Simulation-Program
N-Body simulation program of run-time complexity O(n^2). Different approach to optimize run-time performance are employed.
Program is build upon interface implementation with single level of class hierarchy. Derived classes implement abstract class pure 
virtual functions.
Every class instance has different data container types which are based on STL std::vector, std::valarray and std::array<T,N>, moreover
dynamic allocation by "naked" new and _mm_malloc is in use.
There is an additional set of classes which are based upon AOS data layout which is known to be inferior to SOA 
in terms of vectorisation. The same optimizatio attempts (mainly OpenMP) MT are in place in order to test one solution i.e SOA against
the other i.e AOS.
Extensive profiling is based on Intel VTune Amplifier.
